#!/home/ehsan/anaconda3/bin/python3

import sqlite3
import datetime
import random
import sys
import logging
from config import EMAIL_ACCOUNTS, DB_FILE

def create_tables(conn, table_name):
    conn.execute('PRAGMA foreign_keys = ON;')
    cursor = conn.cursor()

    # Create main table if it doesn't exist
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS "{table_name}" (
            "ID" INTEGER PRIMARY KEY,
            "Professor" TEXT,
            "University" TEXT,
            "Email" TEXT,
            "Website" TEXT
            -- Add other columns as needed
        )
    ''')
    conn.commit()

    # Create email_accounts table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS email_accounts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            from_email TEXT UNIQUE,
            username TEXT,
            password TEXT,
            host TEXT,
            port INTEGER,
            use_ssl BOOLEAN
        )
    ''')
    conn.commit()

    # Insert email accounts from config.py into email_accounts table
    for account in EMAIL_ACCOUNTS:
        cursor.execute('''
            INSERT OR IGNORE INTO email_accounts
            (from_email, username, password, host, port, use_ssl)
            VALUES (?, ?, ?, ?, ?, ?)
        ''', (
            account['from_email'],
            account['username'],
            account['password'],
            account['host'],
            account['port'],
            account['use_ssl']
        ))
    conn.commit()

    # Create search_style_dict table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS search_style_dict (
            search_style INTEGER PRIMARY KEY,
            description TEXT
        )
    ''')
    conn.commit()

    # Insert possible values into search_style_dict
    search_styles = [
        (1, 'breadth-first'),
        (2, 'depth-first'),
        # Add more if needed
    ]
    cursor.executemany('''
        INSERT OR IGNORE INTO search_style_dict (search_style, description)
        VALUES (?, ?)
    ''', search_styles)
    conn.commit()

    # Create answer_state_dict table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS answer_state_dict (
            answer_state INTEGER PRIMARY KEY,
            description TEXT
        )
    ''')
    conn.commit()

    # Insert possible values into answer_state_dict
    answer_states = [
        (0, 'no response yet'),
        (1, 'rejected explicitly'),
        (2, 'no position'),
        (3, 'no position but positive attitude'),
        (4, 'cold answer'),
        (5, 'encouraging answer'),
        (6, 'interview'),
        (10, 'an error occurred'),
    ]
    cursor.executemany('''
        INSERT OR IGNORE INTO answer_state_dict (answer_state, description)
        VALUES (?, ?)
    ''', answer_states)
    conn.commit()

    # Create university_table
    professors_columns = ', '.join([f'"professor{i}" INTEGER' for i in range(1, 11)])
    foreign_keys = ', '.join([f'FOREIGN KEY("professor{i}") REFERENCES "{table_name}"("ID")' for i in range(1, 11)])
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS "university_table" (
            "ID_university" INTEGER PRIMARY KEY AUTOINCREMENT,
            "University" TEXT UNIQUE,
            {professors_columns},
            {foreign_keys}
        )
    ''')
    conn.commit()

    # Create {table_name}_chronology table
    chronology_table = f"{table_name}_chronology"
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS "{chronology_table}" (
            "ID" INTEGER PRIMARY KEY,
            "search_type" INTEGER,
            "search_depth" INTEGER,
            "search_date" INTEGER,
            "data_gathering_completed" BOOLEAN DEFAULT FALSE,
            "data_filtering_completed" BOOLEAN DEFAULT FALSE,
            "html_generation_completed" BOOLEAN DEFAULT FALSE,
            "cv_generation_completed" BOOLEAN DEFAULT FALSE,
            "email_sent" BOOLEAN DEFAULT FALSE,
            "from_email" TEXT,
            "sending_method" INTEGER,
            "answer_state" INTEGER DEFAULT 0,
            "reminder_interval_1" INTEGER,
            "reminder1" INTEGER,
            "reminder_interval_2" INTEGER,
            "reminder2" INTEGER,
            "reminder_interval_3" INTEGER,
            "reminder3" INTEGER,
            FOREIGN KEY("ID") REFERENCES "{table_name}"("ID"),
            FOREIGN KEY("search_type") REFERENCES "search_style_dict"("search_style"),
            FOREIGN KEY("answer_state") REFERENCES "answer_state_dict"("answer_state"),
            FOREIGN KEY("from_email") REFERENCES "email_accounts"("from_email")
        )
    ''')
    conn.commit()

    return chronology_table

def select_professor(cursor, table_name, chronology_table):
    # Enable foreign key constraints
    cursor.execute("PRAGMA foreign_keys = ON;")

    # Get IDs of professors not yet in the chronology table or with answer_state = 0
    cursor.execute(f"""
        SELECT "{table_name}"."ID", "{table_name}"."University"
        FROM "{table_name}"
        LEFT JOIN "{chronology_table}" ON "{table_name}"."ID" = "{chronology_table}"."ID"
        WHERE "{chronology_table}"."ID" IS NULL
    """)
    candidates = cursor.fetchall()
    if not candidates:
        logging.info("No professors available for processing.")
        return None  # Return None to indicate no suitable professor found

    # Shuffle the list to pick a random candidate
    random.shuffle(candidates)

    for candidate in candidates:
        professor_id, university = candidate

        # Check university status
        cursor.execute("""
            SELECT * FROM "university_table" WHERE "University" = ?
        """, (university,))
        university_record = cursor.fetchone()

        if not university_record:
            # University not in university_table, proceed to select professor
            return professor_id, university, 1  # Position 1 (professor1)
        else:
            # Check if another professor can be added
            for i in range(1, 11):
                professor_column = f"professor{i}"
                if university_record[professor_column] is None:
                    # Additional checks can be added here
                    return professor_id, university, i  # Position i (professor2 to professor10)
            # University slots are full, skip
            continue

    logging.info("No suitable professor found.")
    return None  # Return None to indicate no suitable professor found

def add_entries(cursor, table_name, chronology_table, professor_id, university, position):
    # Enable foreign key constraints
    cursor.execute('PRAGMA foreign_keys = ON;')

    # Add entry to chronology table
    today = int(datetime.datetime.now().timestamp())
    cursor.execute(f'''
        INSERT INTO "{chronology_table}" (
            "ID", "search_type", "search_depth", "search_date",
            "data_gathering_completed",
            "data_filtering_completed",
            "html_generation_completed",
            "cv_generation_completed"
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ''', (professor_id, None, None, today, False, False, False, False))

    # Add or update entry in university_table
    cursor.execute('''
        SELECT * FROM "university_table" WHERE "University" = ?
    ''', (university,))
    university_record = cursor.fetchone()

    if not university_record:
        # Insert new university
        columns = '"University", "professor1"'
        values = '?, ?'
        data = (university, professor_id)
        cursor.execute(f'''
            INSERT INTO "university_table" ({columns}) VALUES ({values})
        ''', data)
    else:
        # Update existing university
        professor_column = f"professor{position}"
        cursor.execute(f'''
            UPDATE "university_table"
            SET "{professor_column}" = ?
            WHERE "University" = ?
        ''', (professor_id, university))

    # Commit the changes
    cursor.connection.commit()

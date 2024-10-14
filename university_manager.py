# university_manager.py

import sqlite3
import datetime

def add_professor_to_university(conn, university_name, professor_id, chronology_table):
    cursor = conn.cursor()

    # Check if the university exists
    cursor.execute('''
        SELECT ID_university, professor1, professor2, professor3, professor4, professor5, professor6, professor7, professor8, professor9, professor10
        FROM university_table
        WHERE University = ?
    ''', (university_name,))
    result = cursor.fetchone()

    if result:
        # University exists
        ID_university = result[0]
        professor_columns = result[1:]
        # Find the next available professor slot
        for idx, prof_id in enumerate(professor_columns):
            if prof_id is None:
                column_name = f'professor{idx + 1}'
                cursor.execute(f'''
                    UPDATE university_table
                    SET {column_name} = ?
                    WHERE ID_university = ?
                ''', (professor_id, ID_university))
                conn.commit()
                return True
        # No available slot
        return False
    else:
        # University does not exist, create it and add the professor
        cursor.execute('''
            INSERT INTO university_table (University, professor1)
            VALUES (?, ?)
        ''', (university_name, professor_id))
        conn.commit()
        return True

def can_select_new_professor(conn, university_name, chronology_table):
    cursor = conn.cursor()

    # Check if the university exists
    cursor.execute('''
        SELECT ID_university, professor1, professor2, professor3, professor4, professor5, professor6, professor7, professor8, professor9, professor10
        FROM university_table
        WHERE University = ?
    ''', (university_name,))
    result = cursor.fetchone()

    if not result:
        # University does not exist, so we can select a professor
        return True

    # University exists
    professor_ids = [prof_id for prof_id in result[1:] if prof_id is not None]

    for professor_id in professor_ids:
        # Get the answer_status and send_date for each professor
        cursor.execute(f'''
            SELECT answer_status, send_date
            FROM {chronology_table}
            WHERE ID = ?
        ''', (professor_id,))
        res = cursor.fetchone()
        if res:
            answer_status, send_date = res
            days_since_sent = (datetime.datetime.now() - datetime.datetime.fromtimestamp(send_date)).days if send_date else None

            if answer_status == 0 and days_since_sent is not None and days_since_sent < 7:
                # Professor has not replied, and it's been less than 7 days
                return False
            elif answer_status not in [0, 1, 2, 10]:
                # Professor's answer_status is not in the allowed statuses
                return False
        else:
            # No entry in chronology_table; treat as not allowed
            return False

    # All previous professors meet the conditions
    return True

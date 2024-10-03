#!/home/ehsan/anaconda3/bin/python3

import argparse
import data_gathering
import data_filtering
import modifier
import database_utils
import send_email
import sqlite3
import random
import time
import os
from config import (
    TABLE_NAME,
    SEARCH_DEPTH,
    PROJECT_DIRECTORY,
    DB_FILE,
    TEST_RUN,
    TEST_EMAIL,
    EMAIL_ACCOUNTS,
    REMINDER_INTERVAL_1,
    REMINDER_INTERVAL_2,
    REMINDER_INTERVAL_3,
    SENDING_METHOD,
    SEARCH_STYLE,
    ANSWER_STATE,
    # Add other configurations as needed
)
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Main script to orchestrate the project.')
    parser.add_argument('-i', '--input',
                        help='SQLite database file.')
    parser.add_argument('-t', '--table-name',
                        help='Name of the table in the database.')
    parser.add_argument('-d', '--search-depth', type=int,
                        help='Depth of the link search.')
    parser.add_argument('-p', '--project-directory',
                        help='Project directory for storing data.')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Use command-line arguments if provided; otherwise, use config.py values
    db_file = args.input if args.input else DB_FILE
    table_name = args.table_name if args.table_name else TABLE_NAME
    project_directory = args.project_directory if args.project_directory else PROJECT_DIRECTORY
    search_depth = args.search_depth if args.search_depth else SEARCH_DEPTH

    conn = sqlite3.connect(db_file)
    conn.execute('PRAGMA foreign_keys = ON;')
    cursor = conn.cursor()

    # Create necessary tables
    chronology_table = database_utils.create_tables(conn, table_name)

    # Main loop to process professors
    # Fetch all professors
    cursor.execute(f'''
        SELECT p."ID", p."Professor", p."Email", p."University"
        FROM "{table_name}" p
    ''')

    professors = cursor.fetchall()

    for professor in professors:
        professor_id, professor_name, professor_email, university = professor

        # Check if the professor already exists in the chronology table
        cursor.execute(f'SELECT * FROM "{chronology_table}" WHERE "ID" = ?', (professor_id,))
        chronology_entry = cursor.fetchone()

        if chronology_entry is None:
            # **Determine position for the professor in the university_table**
            position = 1  # Starting with position 1 for simplicity

            # Add entries to chronology_table and university_table
            database_utils.add_entries(cursor, table_name, chronology_table, professor_id, university, position)
            conn.commit()
            cursor.execute(f'SELECT * FROM "{chronology_table}" WHERE "ID" = ?', (professor_id,))
            chronology_entry = cursor.fetchone()

        # Unpack chronology_entry
        (
            c_ID, search_type, c_search_depth, c_search_date,
            data_gathering_completed, data_filtering_completed,
            html_generation_completed, cv_generation_completed,
            email_sent, from_email, sending_method,
            answer_state, reminder_interval_1, reminder1,
            reminder_interval_2, reminder2, reminder_interval_3, reminder3
        ) = chronology_entry

        # Skip if email_sent is True
        if email_sent:
            print(f"Email already sent to Professor ID {professor_id}: {professor_name}")
            continue

        print(f"Processing Professor ID {professor_id}: {professor_name}")

        # Data Gathering
        if not data_gathering_completed:
            # Call data_gathering module with professor_id
            data_gathering.main(db_file, table_name, project_directory, search_depth, professor_id)
            # Update the chronology table
            cursor.execute(f'''
                UPDATE "{chronology_table}"
                SET "data_gathering_completed" = TRUE
                WHERE "ID" = ?
            ''', (professor_id,))
            conn.commit()

        # Data Filtering
        if not data_filtering_completed:
            # Call data_filtering module
            data_filtering.filter_professor_data(professor_name, project_directory)
            # Update the chronology table
            cursor.execute(f'''
                UPDATE "{chronology_table}"
                SET "data_filtering_completed" = TRUE
                WHERE "ID" = ?
            ''', (professor_id,))
            conn.commit()

        # Modifier
        if not html_generation_completed or not cv_generation_completed:
            # Call modifier module
            modifier.modify_template(db_file, table_name, project_directory, professor_id)
            # Update the chronology table
            cursor.execute(f'''
                UPDATE "{chronology_table}"
                SET "html_generation_completed" = TRUE,
                    "cv_generation_completed" = TRUE
                WHERE "ID" = ?
            ''', (professor_id,))
            conn.commit()

        # Email Sending
        # Select a random from_email from EMAIL_ACCOUNTS
        email_account = random.choice(EMAIL_ACCOUNTS)

        from_email = email_account['from_email']

        # Update the chronology table with from_email
        cursor.execute(f'''
            UPDATE "{chronology_table}"
            SET "from_email" = ?
            WHERE "ID" = ?
        ''', (from_email, professor_id))
        conn.commit()

        # Prepare the email content and attachments
        # Load email1.html
        safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
        professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
        email_html_path = os.path.join(professor_dir, 'email1.html')
        cv_pdf_path = os.path.join(professor_dir, 'Ehsan_Ghavimehr_CV_modified.pdf')

        if not os.path.exists(email_html_path):
            print(f"Email HTML file not found for {professor_name}")
            continue
        if not os.path.exists(cv_pdf_path):
            print(f"CV PDF file not found for {professor_name}")
            continue

        with open(email_html_path, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Prepare to send email
        to_email = professor_email
        subject = 'Prospective Ph.D. Student Interested in Cognitive Psychology Research'

        # Handle TEST_RUN
        if TEST_RUN:
            to_email = TEST_EMAIL

        # Send email
        send_email.send_email_smtp(
            to_email=to_email,
            subject=subject,
            html_content=html_content,
            attachment_paths=[cv_pdf_path],
            smtp_config=email_account
        )

        # Update the chronology table
        cursor.execute(f'''
            UPDATE "{chronology_table}"
            SET "email_sent" = TRUE
            WHERE "ID" = ?
        ''', (professor_id,))
        conn.commit()

        # Wait for a random interval between emails
        time_interval = random.randint(60, 180)  # Wait between 1 and 3 minutes
        print(f"Waiting for {time_interval} seconds before sending the next email...")
        time.sleep(time_interval)

    conn.close()

if __name__ == '__main__':
    main()

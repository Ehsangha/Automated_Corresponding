#!/home/ehsan/anaconda3/bin/python3

# main.py

import argparse
import sqlite3
import os
import logging
import sys
import time

# Importing custom modules
import data_gathering
import data_filtering
import modifier
import send_email

# Import configuration settings
from config import (
    TABLE_NAME,
    SEARCH_DEPTH,
    PROJECT_DIRECTORY,
    DB_FILE,
    EMAIL_ACCOUNTS,
    TEST_RUN,
    TEST_EMAIL
)

def setup_logging():
    """
    Configure logging to display messages on the terminal with timestamps and log levels.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

def parse_arguments():
    """
    Parse command-line arguments to allow overriding default configurations.
    """
    parser = argparse.ArgumentParser(
        description='Main script to orchestrate the project workflow.')
    parser.add_argument('-i', '--input',
                        help='SQLite database file.',
                        default=DB_FILE)
    parser.add_argument('-t', '--table-name',
                        help='Name of the table in the database.',
                        default=TABLE_NAME)
    parser.add_argument('-d', '--search-depth', type=int,
                        help='Depth of the link search.',
                        default=SEARCH_DEPTH)
    parser.add_argument('-p', '--project-directory',
                        help='Project directory for storing data.',
                        default=PROJECT_DIRECTORY)
    return parser.parse_args()

def main():
    """
    Main function to orchestrate the entire workflow:
    1. Data Gathering
    2. Data Filtering
    3. Template Modification
    4. Email Sending
    """
    setup_logging()
    args = parse_arguments()

    db_file = args.input
    table_name = args.table_name
    project_directory = args.project_directory
    search_depth = args.search_depth

    # Validate database file
    if not os.path.isfile(db_file):
        logging.error(f"Database file '{db_file}' does not exist.")
        sys.exit(1)

    # Connect to SQLite database
    try:
        conn = sqlite3.connect(db_file)
        conn.row_factory = sqlite3.Row  # Enable row access by column name
        cursor = conn.cursor()
        logging.info(f"Connected to database '{db_file}' successfully.")
    except sqlite3.Error as e:
        logging.error(f"Database connection error: {e}")
        sys.exit(1)

    # Ensure chronology table exists
    chronology_table = f"{table_name}_chronology"
    try:
        cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS "{chronology_table}" (
                "ID" INTEGER PRIMARY KEY,
                "data_gathering_completed" BOOLEAN DEFAULT FALSE,
                "web_scraping_completed" BOOLEAN DEFAULT FALSE,
                "data_filtering_completed" BOOLEAN DEFAULT FALSE,
                "html_generation_completed" BOOLEAN DEFAULT FALSE,
                "cv_generation_completed" BOOLEAN DEFAULT FALSE,
                "email_sent" BOOLEAN DEFAULT FALSE
            );
        """)
        conn.commit()
        logging.info(f"Ensured chronology table '{chronology_table}' exists.")
    except sqlite3.Error as e:
        logging.error(f"Error creating chronology table: {e}")
        conn.close()
        sys.exit(1)

    # Fetch all professors from the main table
    try:
        cursor.execute(f'''
            SELECT "ID", "Professor", "University", "Email"
            FROM "{table_name}"
        ''')
        professors = cursor.fetchall()
        logging.info(f"Fetched {len(professors)} professors from the table '{table_name}'.")
    except sqlite3.Error as e:
        logging.error(f"Error fetching professors: {e}")
        conn.close()
        sys.exit(1)

    if not professors:
        logging.warning("No professors found in the database. Exiting.")
        conn.close()
        sys.exit(0)

    for professor in professors:
        professor_id, professor_name, university, professor_email = professor
        logging.info(f"--- Processing Professor ID {professor_id}: {professor_name} ---")

        # Check chronology table for the professor
        try:
            cursor.execute(f'''
                SELECT 
                    "data_gathering_completed",
                    "web_scraping_completed",
                    "data_filtering_completed",
                    "html_generation_completed",
                    "cv_generation_completed",
                    "email_sent"
                FROM "{chronology_table}"
                WHERE "ID" = ?
            ''', (professor_id,))
            chronology = cursor.fetchone()

            if chronology:
                (data_gathering_completed,
                 web_scraping_completed,
                 data_filtering_completed,
                 html_generation_completed,
                 cv_generation_completed,
                 email_sent) = chronology
            else:
                # If no entry exists, initialize it
                cursor.execute(f'''
                    INSERT INTO "{chronology_table}" (
                        "ID",
                        "data_gathering_completed",
                        "web_scraping_completed",
                        "data_filtering_completed",
                        "html_generation_completed",
                        "cv_generation_completed",
                        "email_sent"
                    ) VALUES (?, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
                ''', (professor_id,))
                conn.commit()
                (data_gathering_completed,
                 web_scraping_completed,
                 data_filtering_completed,
                 html_generation_completed,
                 cv_generation_completed,
                 email_sent) = (False, False, False, False, False, False)
                logging.info(f"Initialized chronology entry for Professor ID {professor_id}.")

        except sqlite3.Error as e:
            logging.error(f"Error accessing chronology table: {e}")
            continue

        # Step 1: Data Gathering
        if not data_gathering_completed:
            logging.info("Starting Data Gathering...")
            try:
                # Call the data_gathering module's main function
                # It returns professor_id and professor_name
                gathered_professor_id, gathered_professor_name = data_gathering.main(
                    db_file=db_file,
                    table_name=table_name,
                    project_directory=project_directory,
                    search_depth=search_depth
                )
                # Update chronology
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "data_gathering_completed" = TRUE
                    WHERE "ID" = ?
                ''', (professor_id,))
                conn.commit()
                logging.info("Data Gathering completed and chronology updated.")
            except Exception as e:
                logging.error(f"Data Gathering failed for Professor ID {professor_id}: {e}")
                continue
        else:
            logging.info("Data Gathering already completed. Skipping...")

        # Step 2: Web Scraping
        if not web_scraping_completed:
            logging.info("Starting Web Scraping...")
            try:
                # Call the web_scraper module's download_webpages function
                # Assuming web_scraper.py is in the same directory and has download_webpages defined
                # data_gathering.main() already called download_webpages, so this step might be redundant
                # If you need to call it separately, uncomment the following lines:
                #
                # import web_scraper
                # web_scraper.download_webpages(cursor, table_name, professor_id, project_directory, search_depth)
                #
                # For now, we'll assume it's handled in data_gathering.main()
                
                # Update chronology
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "web_scraping_completed" = TRUE
                    WHERE "ID" = ?
                ''', (professor_id,))
                conn.commit()
                logging.info("Web Scraping completed and chronology updated.")
            except Exception as e:
                logging.error(f"Web Scraping failed for Professor ID {professor_id}: {e}")
                continue
        else:
            logging.info("Web Scraping already completed. Skipping...")

        # Step 3: Data Filtering
        if not data_filtering_completed:
            logging.info("Starting Data Filtering...")
            try:
                # Call the data_filtering module's filter_professor_data function
                data_filtering.filter_professor_data(professor_name, project_directory)
                # Update chronology
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "data_filtering_completed" = TRUE
                    WHERE "ID" = ?
                ''', (professor_id,))
                conn.commit()
                logging.info("Data Filtering completed and chronology updated.")
            except Exception as e:
                logging.error(f"Data Filtering failed for Professor ID {professor_id}: {e}")
                continue
        else:
            logging.info("Data Filtering already completed. Skipping...")

        # Step 4: Template Modification
        if not (html_generation_completed and cv_generation_completed):
            logging.info("Starting Template Modification...")
            try:
                # Call the modifier module's modify_template function
                modifier.modify_template(
                    db_file=db_file,
                    table_name=table_name,
                    project_directory=project_directory,
                    professor_id=professor_id
                )
                # Update chronology
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "html_generation_completed" = TRUE,
                        "cv_generation_completed" = TRUE
                    WHERE "ID" = ?
                ''', (professor_id,))
                conn.commit()
                logging.info("Template Modification completed and chronology updated.")
            except Exception as e:
                logging.error(f"Template Modification failed for Professor ID {professor_id}: {e}")
                continue
        else:
            logging.info("Template Modification already completed. Skipping...")

        # Step 5: Data Filtering (if needed again)
        # Depending on your workflow, you might not need to run data_filtering again here.

        # Step 6: Email Sending
        if not email_sent:
            logging.info("Starting Email Sending...")
            try:
                # Define paths
                safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
                professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
                email_html_path = os.path.join(professor_dir, 'email1.html')
                cv_pdf_path = os.path.join(professor_dir, 'Ehsan_Ghavimehr_CV_modified.pdf')

                # Verify that the modified email and CV exist
                if not os.path.exists(email_html_path):
                    logging.error(f"Email HTML file not found at {email_html_path}. Skipping email sending.")
                    continue
                if not os.path.exists(cv_pdf_path):
                    logging.error(f"Modified CV PDF file not found at {cv_pdf_path}. Skipping email sending.")
                    continue

                with open(email_html_path, 'r', encoding='utf-8') as f:
                    html_content = f.read()

                # Prepare email details
                to_email = TEST_EMAIL if TEST_RUN else professor_email
                subject = 'Research Collaboration Opportunity'

                attachment_paths = [cv_pdf_path]

                # Select SMTP account (assuming first account)
                smtp_config = EMAIL_ACCOUNTS[0]

                # Send email
                send_email.send_email_smtp(
                    to_email=to_email,
                    subject=subject,
                    html_content=html_content,
                    attachment_paths=attachment_paths,
                    smtp_config=smtp_config
                )

                # Update chronology
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "email_sent" = TRUE
                    WHERE "ID" = ?
                ''', (professor_id,))
                conn.commit()
                logging.info("Email sent and chronology updated.")
            except Exception as e:
                logging.error(f"Email Sending failed for Professor ID {professor_id}: {e}")
                continue
        else:
            logging.info("Email already sent. Skipping...")

        # Optional: Wait between processing professors to avoid rate limits
        time.sleep(2)  # Sleep for 2 seconds

    # Close the database connection
    conn.close()
    logging.info("All professors processed. Database connection closed.")

if __name__ == '__main__':
    main()

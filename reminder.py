# reminder.py

import sqlite3
import datetime
import send_email  # Assuming send_email.py is in the same directory
import logging
import os
import email
from email.utils import formataddr, parsedate_to_datetime
from email import policy
import imaplib
import config  # Import your config.py
import time

def send_reminders(db_file, table_name, project_directory, email_account_id):
    # Configure logging
    logger = logging.getLogger('reminder')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Remove any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler for logging
    fh = logging.FileHandler('reminder.log')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Stream handler for console output
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    chronology_table = f"{table_name}_chronology"

    # Fetch email account details
    cursor.execute('''
        SELECT "from_email", "username", "password", "smtp_host", "smtp_port", "imap_host", "imap_port", "ssl"
        FROM email_accounts
        WHERE "ID" = ?
    ''', (email_account_id,))
    email_account = cursor.fetchone()
    if not email_account:
        logger.error(f"No email account found with ID {email_account_id}")
        conn.close()
        return

    from_email, username, password, smtp_host, smtp_port, imap_host, imap_port, ssl = email_account

    # Fetch professors who need reminders
    professors_to_remind = get_professors_to_remind(cursor, chronology_table, table_name)

    logger.info(f"Found {len(professors_to_remind)} professors to send reminders.")

    for professor in professors_to_remind:
        professor_id, professor_name, professor_email, reminder_number = professor

        logger.info(f"Sending reminder {reminder_number} to Professor ID {professor_id}: {professor_name}")

        # Fetch the original message ID for threading
        message_id, thread_id = fetch_original_message_id(
            imap_host, imap_port, username, password, professor_email, ssl, logger
        )

        if not message_id:
            logger.warning(f"Could not find original message ID for Professor ID {professor_id}. Skipping.")
            continue

        # Generate reminder email content
        email_html_path = os.path.join(
            project_directory, 'data', ''.join(c if c.isalnum() else '_' for c in professor_name), 'email1.html'
        )

        if not os.path.exists(email_html_path):
            logger.warning(f"Email HTML file not found for {professor_name}")
            continue
        else:
            logger.info(f"Email HTML file found for {professor_name}: {email_html_path}")

        with open(email_html_path, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Prepare to send email
        to_email = professor_email
        subject = 'Re: Research Collaboration Opportunity'  # Prefix with 'Re:'

        # Handle TEST_RUN
        if config.TEST_RUN:
            to_email = config.TEST_EMAIL
            logger.info(f"TEST_RUN is enabled. Email will be sent to {config.TEST_EMAIL} instead of {professor_email}")

        # Send reminder email as a reply
        email_sent = send_email.send_email_smtp(
            db_file=db_file,
            email_account_id=email_account_id,
            to_email=to_email,
            subject=subject,
            html_content=html_content,
            attachment_paths=[],
            in_reply_to=message_id,
            references=thread_id
        )

        if email_sent:
            # Update the chronology table
            current_timestamp = int(datetime.datetime.now().timestamp())
            if reminder_number == 1:
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "reminder1" = ?
                    WHERE "ID" = ?
                ''', (current_timestamp, professor_id))
            elif reminder_number == 2:
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "reminder2" = ?
                    WHERE "ID" = ?
                ''', (current_timestamp, professor_id))
            elif reminder_number == 3:
                cursor.execute(f'''
                    UPDATE "{chronology_table}"
                    SET "reminder3" = ?
                    WHERE "ID" = ?
                ''', (current_timestamp, professor_id))
            conn.commit()
            logger.info(f"Reminder {reminder_number} sent to {to_email} for Professor ID {professor_id}")
        else:
            logger.error(f"Failed to send reminder {reminder_number} to {to_email} for Professor ID {professor_id}")

    conn.close()

def get_professors_to_remind(cursor, chronology_table, table_name):
    professors = []
    current_timestamp = int(datetime.datetime.now().timestamp())

    # Fetch professors who need reminder 1
    cursor.execute(f'''
        SELECT c."ID", p."Professor", p."Email"
        FROM "{chronology_table}" c
        JOIN "{table_name}" p ON c."ID" = p."ID"
        WHERE c."email_sent" = 1
          AND c."reminder1" IS NULL
          AND c."answer_status" IS NULL
          AND (? - c."send_date") >= ?
    ''', (current_timestamp, config.REMINDER_INTERVAL_1 * 86400))
    professors.extend([(row[0], row[1], row[2], 1) for row in cursor.fetchall()])

    # Fetch professors who need reminder 2
    cursor.execute(f'''
        SELECT c."ID", p."Professor", p."Email"
        FROM "{chronology_table}" c
        JOIN "{table_name}" p ON c."ID" = p."ID"
        WHERE c."reminder1" IS NOT NULL
          AND c."reminder2" IS NULL
          AND c."answer_status" IS NULL
          AND (? - c."reminder1") >= ?
    ''', (current_timestamp, config.REMINDER_INTERVAL_2 * 86400))
    professors.extend([(row[0], row[1], row[2], 2) for row in cursor.fetchall()])

    # Fetch professors who need reminder 3
    cursor.execute(f'''
        SELECT c."ID", p."Professor", p."Email"
        FROM "{chronology_table}" c
        JOIN "{table_name}" p ON c."ID" = p."ID"
        WHERE c."reminder2" IS NOT NULL
          AND c."reminder3" IS NULL
          AND c."answer_status" IS NULL
          AND (? - c."reminder2") >= ?
    ''', (current_timestamp, config.REMINDER_INTERVAL_3 * 86400))
    professors.extend([(row[0], row[1], row[2], 3) for row in cursor.fetchall()])

    return professors

def fetch_original_message_id(imap_host, imap_port, username, password, to_email, ssl, logger):
    try:
        if ssl:
            imap = imaplib.IMAP4_SSL(imap_host, imap_port)
        else:
            imap = imaplib.IMAP4(imap_host, imap_port)
        imap.login(username, password)

        # Select the 'Sent' folder
        if 'gmail' in imap_host.lower():
            folder_name = '"[Gmail]/Sent Mail"'
        elif 'hostinger' in imap_host.lower():
            folder_name = 'INBOX.Sent'
        else:
            folder_name = 'Sent'

        imap.select(folder_name)

        # Search for the sent email to the professor
        status, data = imap.search(None, f'(TO "{to_email}")')
        if status != 'OK':
            logger.warning(f"IMAP search failed for {to_email}")
            imap.logout()
            return None, None

        email_ids = data[0].split()
        if not email_ids:
            logger.warning(f"No sent emails found to {to_email}")
            imap.logout()
            return None, None

        # Fetch the latest email
        latest_email_id = email_ids[-1]
        status, msg_data = imap.fetch(latest_email_id, '(RFC822)')
        if status != 'OK':
            logger.warning(f"Failed to fetch email with ID {latest_email_id}")
            imap.logout()
            return None, None

        raw_email = msg_data[0][1]
        email_message = email.message_from_bytes(raw_email, policy=policy.default)

        message_id = email_message.get('Message-ID')
        references = email_message.get('References') or message_id

        imap.logout()

        return message_id, references

    except Exception as e:
        logger.exception(f"Error fetching original message ID: {e}")
        return None, None

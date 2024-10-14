# send_email.py

import os
import sqlite3
import logging
import html2text
import time
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
import smtplib
import imaplib
import email
import ssl

# Import the gmail library if available
try:
    from gmail import GMail, Message
except ImportError:
    GMail = None  # Handle the case where gmail library is not installed

def send_email_smtp(db_file, email_account_id, to_email, subject, html_content, attachment_paths, in_reply_to=None, references=None):
    # Configure logging
    logger = logging.getLogger('send_email')
    logger.setLevel(logging.DEBUG)  # Set to DEBUG level
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Remove any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler for logging
    fh = logging.FileHandler('send_email.log')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Stream handler for console output
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Connect to the database to retrieve email account settings
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Fetch the email account details
    cursor.execute('''
        SELECT "from_email", "username", "password", "smtp_host", "smtp_port", "imap_host", "imap_port", "ssl"
        FROM email_accounts
        WHERE "ID" = ?
    ''', (email_account_id,))
    result = cursor.fetchone()
    conn.close()

    if not result:
        logger.error(f"No email account found with ID {email_account_id}")
        return False

    from_email, username, app_password, smtp_host, smtp_port, imap_host, imap_port, ssl_flag = result

    # Detect email provider based on from_email domain
    if from_email.lower().endswith('@gmail.com'):
        logger.info(f"Detected Gmail account: {from_email}")
        # Use send_email_gmail function
        return send_email_gmail(
            from_email, username, app_password, to_email, subject, html_content, attachment_paths, in_reply_to, references, logger
        )
    else:
        logger.info(f"Using SMTP to send email from {from_email}")
        # Use send_email_hostinger function
        return send_email_hostinger(
            from_email, username, app_password, smtp_host, smtp_port, imap_host, imap_port, ssl_flag,
            to_email, subject, html_content, attachment_paths, in_reply_to, references, logger
        )

def send_email_gmail(from_email, username, password, to_email, subject, html_content, attachment_paths, in_reply_to=None, references=None, logger=None):
    if GMail is None:
        if logger:
            logger.error("gmail module is not installed. Cannot send email via Gmail.")
        else:
            print("gmail module is not installed. Cannot send email via Gmail.")
        return False

    try:
        # Initialize the Gmail object
        gmail = GMail(username, password)

        # Create a plain text version of the HTML content
        plain_text_content = html2text.html2text(html_content)

        # Create the email message using MIMEMultipart for proper formatting
        msg = MIMEMultipart('mixed')
        msg['Subject'] = subject
        msg['From'] = from_email
        msg['To'] = to_email

        # Add In-Reply-To and References headers if provided
        if in_reply_to:
            msg['In-Reply-To'] = in_reply_to
        if references:
            msg['References'] = references

        # Attach plain text and HTML content
        alternative_part = MIMEMultipart('alternative')
        part1 = MIMEText(plain_text_content, 'plain')
        part2 = MIMEText(html_content, 'html')
        alternative_part.attach(part1)
        alternative_part.attach(part2)
        msg.attach(alternative_part)

        # Attach files
        for path in attachment_paths:
            if os.path.exists(path):
                with open(path, 'rb') as f:
                    part = MIMEApplication(f.read(), Name=os.path.basename(path))
                part['Content-Disposition'] = f'attachment; filename="{os.path.basename(path)}"'
                msg.attach(part)
            else:
                if logger:
                    logger.warning(f"Attachment file not found: {path}")
                else:
                    print(f"Attachment file not found: {path}")

        # Send the email
        gmail.send(msg)
        if logger:
            logger.info(f"Email sent to {to_email} from {from_email}")
        else:
            print(f"Email sent to {to_email} from {from_email}")

        return True

    except Exception as e:
        if logger:
            logger.exception(f"Failed to send email to {to_email}: {e}")
        else:
            print(f"Failed to send email to {to_email}: {e}")
        return False

def send_email_hostinger(from_email, username, password, smtp_host, smtp_port, imap_host, imap_port, ssl_flag, to_email, subject, html_content, attachment_paths, in_reply_to=None, references=None, logger=None):
    try:
        # Create the email message
        msg = MIMEMultipart('mixed')
        msg['From'] = from_email
        msg['To'] = to_email
        msg['Subject'] = subject

        # Add In-Reply-To and References headers if provided
        if in_reply_to:
            msg['In-Reply-To'] = in_reply_to
        if references:
            msg['References'] = references

        # Create a plain text version of the HTML content
        plain_text_content = html2text.html2text(html_content)

        # Attach the plain text and HTML content
        alternative_part = MIMEMultipart('alternative')
        alternative_part.attach(MIMEText(plain_text_content, 'plain'))
        alternative_part.attach(MIMEText(html_content, 'html'))
        msg.attach(alternative_part)

        # Attach files
        for path in attachment_paths:
            if os.path.exists(path):
                with open(path, 'rb') as f:
                    part = MIMEApplication(f.read(), Name=os.path.basename(path))
                part['Content-Disposition'] = f'attachment; filename="{os.path.basename(path)}"'
                msg.attach(part)
            else:
                if logger:
                    logger.warning(f"Attachment file not found: {path}")
                else:
                    print(f"Attachment file not found: {path}")

        # Send the email using SMTP
        context = ssl.create_default_context()
        if ssl_flag:
            server = smtplib.SMTP_SSL(smtp_host, smtp_port, context=context)
        else:
            server = smtplib.SMTP(smtp_host, smtp_port)
            server.starttls(context=context)
        server.login(username, password)
        server.sendmail(from_email, to_email, msg.as_string())
        server.quit()

        if logger:
            logger.info(f"Email sent to {to_email} from {from_email}")
        else:
            print(f"Email sent to {to_email} from {from_email}")

        # Save the sent email to the 'Sent' folder via IMAP
        save_email_to_sent_folder(imap_host, imap_port, username, password, ssl_flag, msg, logger)

        return True

    except Exception as e:
        if logger:
            logger.exception(f"Failed to send email to {to_email}: {e}")
        else:
            print(f"Failed to send email to {to_email}: {e}")
        return False

def save_email_to_sent_folder(imap_host, imap_port, username, password, ssl_flag, msg, logger=None):
    try:
        if ssl_flag:
            imap = imaplib.IMAP4_SSL(imap_host, imap_port)
        else:
            imap = imaplib.IMAP4(imap_host, imap_port)
        imap.login(username, password)

        if 'gmail' in imap_host.lower():
            folder_name = '"[Gmail]/Sent Mail"'
        elif 'hostinger' in imap_host.lower():
            folder_name = 'INBOX.Sent'
        else:
            folder_name = 'Sent'

        imap.select(folder_name)
        imap.append(folder_name, '\\Seen', imaplib.Time2Internaldate(time.time()), msg.as_bytes())
        imap.logout()

        if logger:
            logger.info("Email saved to Sent folder")
        else:
            print("Email saved to Sent folder")

    except Exception as e:
        if logger:
            logger.exception(f"Failed to save email to Sent folder: {e}")
        else:
            print(f"Failed to save email to Sent folder: {e}")

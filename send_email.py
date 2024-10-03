# send_email.py

import smtplib
import ssl
import imaplib
from email.message import EmailMessage
import logging
import os
from config import EMAIL_ACCOUNTS

def send_email_smtp(to_email, subject, html_content, attachment_paths, smtp_config):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = smtp_config['from_email']
    msg['To'] = to_email
    msg.set_content('This is an HTML email. If you see this, your email client does not support HTML.')
    msg.add_alternative(html_content, subtype='html')

    # Attach files
    for path in attachment_paths:
        if not os.path.exists(path):
            logging.warning(f"Attachment {path} does not exist.")
            continue
        with open(path, 'rb') as f:
            file_data = f.read()
            file_name = os.path.basename(path)
            maintype, subtype = 'application', 'octet-stream'
            msg.add_attachment(file_data, maintype=maintype, subtype=subtype, filename=file_name)

    # Connect to SMTP server and send email
    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL(smtp_config['host'], smtp_config['port'], context=context) as server:
            server.login(smtp_config['username'], smtp_config['password'])
            server.send_message(msg)
        logging.info(f"Email sent to {to_email} successfully.")
    except Exception as e:
        logging.error(f"Failed to send email to {to_email}: {e}")
        return

    # Save the sent email to the "Sent" folder via IMAP
    save_to_sent_imap(msg, smtp_config)

def save_to_sent_imap(email_message, smtp_config):
    try:
        # IMAP Configuration
        imap_host = smtp_config.get('imap_host', 'imap.hostinger.com')  # Default Hostinger IMAP
        imap_port = smtp_config.get('imap_port', 993)
        imap_user = smtp_config['username']
        imap_pass = smtp_config['password']

        # Connect to IMAP server
        with imaplib.IMAP4_SSL(imap_host, imap_port) as imap:
            imap.login(imap_user, imap_pass)
            imap.select('"[Gmail]/Sent Mail"')  # Modify if Hostinger uses a different "Sent" folder name

            # Prepare the email for appending
            raw_email = email_message.as_bytes()

            # Append the email to the "Sent" folder
            imap.append('"Sent"', '\\Seen', imaplib.Time2Internaldate(time.time()), raw_email)
            logging.info("Email appended to 'Sent' folder successfully.")
    except Exception as e:
        logging.error(f"Failed to append email to 'Sent' folder: {e}")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Send an email and save it to the Sent folder.')
    parser.add_argument('-to', '--to_email', required=True, help='Recipient email address.')
    parser.add_argument('-s', '--subject', required=True, help='Subject of the email.')
    parser.add_argument('-c', '--content', required=True, help='HTML content of the email.')
    parser.add_argument('-a', '--attachments', nargs='*', default=[], help='Paths to attachment files.')
    parser.add_argument('-smtp', '--smtp_account', type=int, default=0, help='Index of the SMTP account to use from config.py')
    args = parser.parse_args()

    # Select SMTP account based on the provided index
    try:
        smtp_config = EMAIL_ACCOUNTS[args.smtp_account]
    except IndexError:
        logging.error(f"SMTP account index {args.smtp_account} is out of range.")
        sys.exit(1)

    send_email_smtp(args.to_email, args.subject, args.content, args.attachments, smtp_config)

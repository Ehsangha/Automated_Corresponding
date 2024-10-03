#!/home/ehsan/anaconda3/bin/python3

# web_scraper.py

import os
import requests
from urllib.parse import urljoin, urlparse
from bs4 import BeautifulSoup
import datetime
import logging
import sqlite3
from config import DB_FILE, TABLE_NAME, PROJECT_DIRECTORY, SEARCH_DEPTH

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def download_webpages(cursor, table_name, professor_id, project_directory, search_depth=1):
    # Get professor details
    cursor.execute(f"""
        SELECT "Professor", "Website" FROM "{table_name}" WHERE "ID" = ?
    """, (professor_id,))
    result = cursor.fetchone()
    if not result:
        logging.error(f"No professor found with ID {professor_id}")
        return
    professor_name, website_url = result

    # Create directory
    safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
    professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
    os.makedirs(professor_dir, exist_ok=True)

    # Start downloading
    visited_urls = set()
    to_visit = [(website_url, 0)]
    allowed_tlds = {'edu', 'org', 'ac', 'net', 'gov', 'uk', 'ca', 'au'}  # Add other TLDs as needed

    def download_page(current_url, depth):
        if current_url in visited_urls or depth > search_depth:
            return
        visited_urls.add(current_url)
        try:
            response = requests.get(current_url, timeout=10)
            if response.status_code == 200:
                # Save HTML content
                parsed_url = urlparse(current_url)
                file_ext = os.path.splitext(parsed_url.path)[1].lower()
                if file_ext == '.pdf':
                    # Download PDF
                    filename = os.path.join(professor_dir, f"file_{len(visited_urls)}.pdf")
                    with open(filename, 'wb') as f:
                        f.write(response.content)
                else:
                    filename = os.path.join(professor_dir, f"page_{len(visited_urls)}.html")
                    with open(filename, 'w', encoding='utf-8') as f:
                        f.write(response.text)
                    if depth < search_depth:
                        # Parse links
                        soup = BeautifulSoup(response.text, 'html.parser')
                        for link_tag in soup.find_all('a', href=True):
                            link = link_tag['href']
                            full_link = urljoin(current_url, link)
                            parsed_link = urlparse(full_link)
                            link_tld = parsed_link.netloc.split('.')[-1]
                            file_ext = os.path.splitext(parsed_link.path)[1].lower()
                            if file_ext == '.pdf':
                                to_visit.append((full_link, depth + 1))
                            elif link_tld in allowed_tlds and full_link not in visited_urls:
                                to_visit.append((full_link, depth + 1))
            else:
                logging.warning(f"Failed to download {current_url}: HTTP {response.status_code}")
        except requests.RequestException as e:
            logging.error(f"Error downloading {current_url}: {e}")

    while to_visit:
        current_url, depth = to_visit.pop(0)
        download_page(current_url, depth)

    # Update chronology table
    today = (datetime.date.today() - datetime.date(1982, 4, 17)).days
    cursor.execute(f"""
        UPDATE "{table_name}_chronology"
        SET "search_type" = ?, "search_depth" = ?, "search_date" = ?
        WHERE "ID" = ?
    """, (1, search_depth, today, professor_id))
    cursor.connection.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Web Scraper Script')
    parser.add_argument('-i', '--input', help='SQLite database file.')
    parser.add_argument('-t', '--table-name', help='Name of the table in the database.')
    parser.add_argument('-d', '--search-depth', type=int, help='Depth of the link search.')
    parser.add_argument('-p', '--project-directory', help='Project directory for storing data.')
    parser.add_argument('-pid', '--professor-id', type=int, required=True, help='Professor ID to scrape for.')
    args = parser.parse_args()

    db_file = args.input if args.input else DB_FILE
    table_name = args.table_name if args.table_name else TABLE_NAME
    project_directory = args.project_directory if args.project_directory else PROJECT_DIRECTORY
    search_depth = args.search_depth if args.search_depth else SEARCH_DEPTH
    professor_id = args.professor_id

    conn = sqlite3.connect(db_file)
    conn.execute('PRAGMA foreign_keys = ON;')
    cursor = conn.cursor()

    download_webpages(cursor, table_name, professor_id, project_directory, search_depth)

    conn.close()

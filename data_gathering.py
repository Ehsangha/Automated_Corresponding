#!/home/ehsan/anaconda3/bin/python3

import sqlite3
import os
import json
import logging
from scholarly import scholarly
from serpapi import GoogleSearch
from Bio import Entrez
from linkedin_api import Linkedin
from habanero import Crossref
import requests
from datetime import datetime
from config import (
    SERPAPI_API_KEY,
    LINKEDIN_EMAIL,
    LINKEDIN_PASSWORD,
    ORCID_CLIENT_ID,
    ORCID_CLIENT_SECRET,
    DB_FILE,
    TABLE_NAME,
    PROJECT_DIRECTORY,
    SEARCH_DEPTH
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize Crossref
cr = Crossref()

def main(db_file, table_name, project_directory, search_depth, professor_id):
    conn = sqlite3.connect(db_file)
    conn.execute('PRAGMA foreign_keys = ON;')
    cursor = conn.cursor()

    # Fetch professor details using professor_id
    cursor.execute(f'''
        SELECT "Professor", "University", "Email", "Website"
        FROM "{table_name}"
        WHERE "ID" = ?
    ''', (professor_id,))
    result = cursor.fetchone()
    if not result:
        logging.error(f"No professor found with ID {professor_id}")
        conn.close()
        return

    professor_name, university, email, website = result

    # Log the fetched data
    logging.info(f"Professor ID: {professor_id}")
    logging.info(f"Professor Name: {professor_name}")
    logging.info(f"University: {university}")
    logging.info(f"Email: {email}")
    logging.info(f"Website: {website}")

    # Create a safe directory name for the professor
    safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
    professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
    os.makedirs(professor_dir, exist_ok=True)

    # Gather data from various sources
    professor_data = gather_professor_data(professor_name)

    # Save the gathered data to a JSON file
    data_file = os.path.join(professor_dir, 'professor_data.json')
    with open(data_file, 'w', encoding='utf-8') as f:
        json.dump(professor_data, f, ensure_ascii=False, indent=4)
    logging.info(f"Professor data saved to {data_file}")

    # Close the database connection
    conn.close()

def gather_professor_data(professor_name):
    data = {}

    # Google Scholar
    data['Google Scholar'] = search_google_scholar(professor_name)

    # SerpApi Google Scholar
    data['SerpApi Google Scholar'] = search_serpapi_google_scholar(professor_name)

    # PubMed
    data['PubMed'] = search_pubmed(professor_name)

    # LinkedIn
    data['LinkedIn'] = search_linkedin(professor_name)

    # CrossRef
    data['CrossRef'] = search_crossref(professor_name)

    # Semantic Scholar
    data['Semantic Scholar'] = search_semantic_scholar(professor_name)

    # ORCID
    data['ORCID'] = search_orcid(professor_name)

    return data

def search_google_scholar(professor_name):
    try:
        search_query = scholarly.search_author(professor_name)
        author = next(search_query, None)
        if author:
            author = scholarly.fill(author)  # Fetch publications
            publications = []
            for pub in author.get('publications', []):
                pub_filled = scholarly.fill(pub)
                publications.append({
                    'title': pub_filled.get('bib', {}).get('title'),
                    'year': pub_filled.get('bib', {}).get('pub_year'),
                    'abstract': pub_filled.get('bib', {}).get('abstract'),
                    'url': pub_filled.get('bib', {}).get('url'),
                    'citation_count': pub_filled.get('num_citations')
                })
            return publications
        else:
            logging.info(f"No Google Scholar profile found for {professor_name}")
    except Exception as e:
        logging.error(f"Google Scholar search error: {e}")
    return []

def search_serpapi_google_scholar(professor_name):
    params = {
        "engine": "google_scholar",
        "q": professor_name,
        "api_key": SERPAPI_API_KEY
    }
    try:
        search = GoogleSearch(params)
        results = search.get_dict()
        publications = []
        for result in results.get("organic_results", []):
            publications.append({
                'title': result.get('title'),
                'snippet': result.get('snippet'),
                'publication_year': result.get('publication_year'),
                'link': result.get('link'),
                'citations': result.get('citations', {}).get('total')
            })
        return publications
    except Exception as e:
        logging.error(f"SerpAPI error: {e}")
    return []

def search_pubmed(professor_name):
    Entrez.email = "ehsanghavimehr@gmail.com"  # Required by NCBI
    try:
        handle = Entrez.esearch(db="pubmed", term=professor_name, retmax=10)
        record = Entrez.read(handle)
        pmids = record['IdList']
        if pmids:
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
            abstracts = handle.read()
            return abstracts
        else:
            logging.info(f"No PubMed results for {professor_name}")
    except Exception as e:
        logging.error(f"PubMed search error: {e}")
    return ""

def search_linkedin(professor_name):
    try:
        api = Linkedin(LINKEDIN_EMAIL, LINKEDIN_PASSWORD)
        people = api.search_people(keywords=professor_name, limit=1)
        if people:
            profile = people[0]
            return {
                'name': profile.get('public_id'),
                'occupation': profile.get('occupation'),
                'location': profile.get('locationName'),
                'profile_url': f"https://www.linkedin.com/in/{profile.get('public_id')}"
            }
        else:
            logging.info(f"No LinkedIn profile found for {professor_name}")
    except Exception as e:
        logging.error(f"LinkedIn search error: {e}")
    return {}

def search_crossref(professor_name):
    try:
        works = cr.works(query=professor_name, limit=10)
        items = works['message']['items']
        publications = []
        for item in items:
            publications.append({
                'title': item.get('title')[0] if item.get('title') else '',
                'author': [author.get('given') + ' ' + author.get('family') for author in item.get('author', [])],
                'year': item.get('published-print', {}).get('date-parts', [[None]])[0][0],
                'journal': item.get('container-title')[0] if item.get('container-title') else '',
                'doi': item.get('DOI'),
                'URL': item.get('URL')
            })
        return publications
    except Exception as e:
        logging.error(f"CrossRef search error: {e}")
    return []

def search_semantic_scholar(professor_name):
    try:
        url = f"https://api.semanticscholar.org/graph/v1/author/search?query={professor_name}&limit=1"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if data['data']:
                author_id = data['data'][0]['authorId']
                papers_url = f"https://api.semanticscholar.org/graph/v1/author/{author_id}/papers?limit=10&fields=title,year,abstract,externalIds"
                papers_response = requests.get(papers_url)
                if papers_response.status_code == 200:
                    papers_data = papers_response.json()
                    papers = []
                    for paper in papers_data.get('data', []):
                        papers.append({
                            'title': paper.get('title'),
                            'year': paper.get('year'),
                            'abstract': paper.get('abstract'),
                            'doi': paper.get('externalIds', {}).get('DOI')
                        })
                    return papers
        logging.info(f"No Semantic Scholar profile found for {professor_name}")
    except Exception as e:
        logging.error(f"Semantic Scholar search error: {e}")
    return []

def search_orcid(professor_name):
    try:
        url = f"https://pub.orcid.org/v3.0/search/?q={professor_name}"
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            data = response.json()
            if data.get('result'):
                return data['result']
        logging.info(f"No ORCID profile found for {professor_name}")
    except Exception as e:
        logging.error(f"ORCID search error: {e}")
    return []

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Data Gathering Script')
    parser.add_argument('-i', '--input', help='SQLite database file.')
    parser.add_argument('-t', '--table-name', help='Name of the table in the database.')
    parser.add_argument('-d', '--search-depth', type=int, help='Depth of the link search.')
    parser.add_argument('-p', '--project-directory', help='Project directory for storing data.')
    parser.add_argument('-pid', '--professor-id', type=int, required=True, help='Professor ID to gather data for.')
    args = parser.parse_args()

    db_file = args.input if args.input else DB_FILE
    table_name = args.table_name if args.table_name else TABLE_NAME
    project_directory = args.project_directory if args.project_directory else PROJECT_DIRECTORY
    search_depth = args.search_depth if args.search_depth else SEARCH_DEPTH
    professor_id = args.professor_id

    main(db_file, table_name, project_directory, search_depth, professor_id)

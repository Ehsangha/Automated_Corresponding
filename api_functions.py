#!/home/ehsan/anaconda3/bin/python3

# api_functions.py

import logging
import threading
import requests
from scholarly import scholarly
from serpapi import GoogleSearch
from Bio import Entrez
from habanero import Crossref
from config import SERPAPI_KEY, ENTREZ_EMAIL

Entrez.email = ENTREZ_EMAIL
cr = Crossref()

def search_google_scholar(name):
    try:
        search_query = scholarly.search_author(name)
        author = next(search_query, None)
        if author:
            # Fetch author data with only basic sections to reduce data size
            author = scholarly.fill(author, sections=['basics', 'publications'])
            articles = []
            publications = author.get('publications', [])
            # Sort publications by year (if available)
            publications.sort(key=lambda x: int(x.get('bib', {}).get('pub_year', 0)), reverse=True)
            # Process up to 5 recent publications
            for pub in publications[:5]:
                pub_bib = pub.get('bib', {})
                year = pub_bib.get('pub_year')
                if year and int(year) >= 2022:
                    first_author = pub_bib.get('author', '').split(' and ')[0]
                    if name.lower() in first_author.lower():
                        article_data = {
                            'title': pub_bib.get('title'),
                            'abstract': pub_bib.get('abstract', ''),
                            'year': year,
                            'author_position': 'First Author',
                            'pdf_link': pub.get('eprint_url'),
                        }
                        articles.append(article_data)
                else:
                    # Skip older publications
                    continue
            return articles
        else:
            logging.info(f"No Google Scholar profile found for {name}")
    except Exception as e:
        logging.error(f"Google Scholar search error: {e}")
    return []


def search_serpapi_google_scholar(name):
    params = {
        "engine": "google_scholar",
        "q": name,
        "api_key": SERPAPI_KEY  # Use API key from config.py
    }
    try:
        search = GoogleSearch(params)
        results = search.get_dict()
        articles = []
        for result in results.get("organic_results", []):
            year = result.get('publication_date', '')[:4]
            if year and int(year) >= 2022:
                # Assuming the first author is the professor
                article_data = {
                    'title': result.get('title'),
                    'abstract': result.get('snippet'),
                    'year': year,
                    'author_position': 'First Author',
                    'pdf_link': result.get('link'),
                }
                articles.append(article_data)
        return articles
    except Exception as e:
        logging.error(f"SerpAPI error: {e}")
    return []

def search_pubmed(name):
    try:
        # Search for articles by the author in the last 2 years
        current_year = 2024
        start_year = current_year - 2
        search_term = f"{name}[Author] AND ({start_year}:{current_year}[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=20)
        record = Entrez.read(handle)
        pmids = record['IdList']
        articles = []
        if pmids:
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="xml")
            records = Entrez.read(handle)
            for article in records['PubmedArticle']:
                article_data = {}
                article_data['title'] = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_list = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                article_data['abstract'] = ' '.join(abstract_list)
                article_data['year'] = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year')
                # Check if the professor is the first author
                authors = article['MedlineCitation']['Article'].get('AuthorList', [])
                if authors and 'LastName' in authors[0] and name.split()[-1].lower() == authors[0]['LastName'].lower():
                    article_data['author_position'] = 'First Author'
                    # PubMed does not provide PDF links directly
                    article_data['pdf_link'] = None
                    articles.append(article_data)
        else:
            logging.info(f"No PubMed results for {name}")
        return articles
    except Exception as e:
        logging.error(f"PubMed search error: {e}")
    return []

def search_crossref(name):
    try:
        works = cr.works(query=name, limit=10, filter={'from-pub-date': '2022-01-01'})
        articles = []
        for item in works['message']['items']:
            year = item.get('published-print', {}).get('date-parts', [[None]])[0][0]
            if year and int(year) >= 2022:
                first_author = item.get('author', [{}])[0].get('family', '')
                if name.split()[-1].lower() == first_author.lower():
                    article_data = {
                        'title': item.get('title', [])[0] if item.get('title') else None,
                        'abstract': item.get('abstract'),
                        'year': str(year),
                        'author_position': 'First Author',
                        'pdf_link': item.get('URL'),
                    }
                    articles.append(article_data)
        return articles
    except Exception as e:
        logging.error(f"CrossRef search error: {e}")
    return []

def search_orcid(name):
    try:
        url = f"https://pub.orcid.org/v3.0/search/?q=given-names:{name}"
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            data = response.json()
            results = data.get('result', [])
            exact_matches = []
            for result in results:
                orcid_id = result['orcid-identifier']['path']
                personal_details_url = f"https://pub.orcid.org/v3.0/{orcid_id}/personal-details"
                personal_response = requests.get(personal_details_url, headers=headers)
                if personal_response.status_code == 200:
                    personal_data = personal_response.json()
                    given_names = personal_data.get('name', {}).get('given-names', {}).get('value', '').lower()
                    family_name = personal_data.get('name', {}).get('family-name', {}).get('value', '').lower()
                    full_name = f"{given_names} {family_name}"
                    if full_name.strip().lower() == name.strip().lower():
                        # Exact match found
                        # Fetch works
                        works_url = f"https://pub.orcid.org/v3.0/{orcid_id}/works"
                        works_response = requests.get(works_url, headers=headers)
                        if works_response.status_code == 200:
                            works_data = works_response.json()
                            articles = []
                            for group in works_data.get('group', []):
                                summary = group.get('work-summary', [])[0]
                                pub_year = summary.get('publication-date', {}).get('year', {}).get('value')
                                if pub_year and int(pub_year) >= 2022:
                                    article_data = {
                                        'title': summary.get('title', {}).get('title', {}).get('value'),
                                        'abstract': None,  # ORCID doesn't provide abstracts
                                        'year': pub_year,
                                        'author_position': 'Author',  # ORCID doesn't specify author position
                                        'pdf_link': None,  # No direct PDF link
                                    }
                                    articles.append(article_data)
                            return articles
            logging.info(f"No exact ORCID profile found for {name}")
        else:
            logging.info(f"No ORCID profile found for {name}")
    except Exception as e:
        logging.error(f"ORCID search error: {e}")
    return []

def gather_professor_data(professor_name):
    # Gather data from all sources
    data_sources = {
        "Google Scholar": search_google_scholar,
        "SerpAPI": search_serpapi_google_scholar,
        "PubMed": search_pubmed,
        "CrossRef": search_crossref,
        "ORCID": search_orcid
    }

    all_data = {}
    threads = []

    def gather_data(source_name, func):
        data = func(professor_name)
        all_data[source_name] = data

    for source_name, func in data_sources.items():
        thread = threading.Thread(target=gather_data, args=(source_name, func))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    return all_data

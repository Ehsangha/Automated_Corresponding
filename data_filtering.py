#!/home/ehsan/anaconda3/bin/python3

# data_filtering.py

import os
import json
import logging
from datetime import datetime

def filter_professor_data(professor_name, project_directory):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
    data_dir = os.path.join(project_directory, 'data', safe_professor_name)
    data_file = os.path.join(data_dir, 'professor_data.json')
    filtered_file = os.path.join(data_dir, 'professor_data_filtered.json')
    
    if not os.path.exists(data_file):
        logging.error(f"{data_file} does not exist.")
        return
    
    with open(data_file, 'r', encoding='utf-8') as f:
        professor_data = json.load(f)
    
    filtered_articles = []
    current_year = datetime.now().year
    start_year = current_year - 2
    
    for source, articles in professor_data.items():
        for article in articles:
            if isinstance(article, dict):
                # Proceed if article is a dictionary
                article_year = article.get('year')
                if article_year and start_year <= int(article_year) <= current_year:
                    author_position = article.get('author_position')
                    # Assuming 'author_position' field exists; adjust as needed
                    if author_position in ['First Author', 'Corresponding Author']:
                        filtered_articles.append({
                            'title': article.get('title'),
                            'abstract': article.get('abstract'),
                            'year': article_year,
                            'author_position': author_position,
                            'pdf_link': article.get('pdf_link'),
                            'source': source
                        })
            else:
                logging.warning(f"Article from source '{source}' is not in expected format.")
    
    # Sort articles from newest to oldest
    filtered_articles.sort(key=lambda x: int(x['year']), reverse=True)
    
    with open(filtered_file, 'w', encoding='utf-8') as f:
        json.dump(filtered_articles, f, ensure_ascii=False, indent=4)
    
    logging.info(f"Filtered data saved to {filtered_file}")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Data Filtering Script')
    parser.add_argument('-n', '--name', required=True, help='Professor Name')
    parser.add_argument('-p', '--project-directory', help='Project directory for storing data.')
    args = parser.parse_args()

    professor_name = args.name
    project_directory = args.project_directory if args.project_directory else '/home/ehsan/usr/mailling_list'  # Update as needed

    filter_professor_data(professor_name, project_directory)

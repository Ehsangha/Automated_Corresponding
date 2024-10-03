#!/home/ehsan/anaconda3/bin/python3

# modifier.py

import os
import sqlite3
import sys
import json
import subprocess
from openai import OpenAI
from config import (
    OPENAI_API_KEY,
    OPENAI_SUMMARY_MAX_TOKENS,
    OPENAI_SUMMARY_TEMPERATURE,
    OPENAI_PARAGRAPH_MAX_TOKENS,
    OPENAI_PARAGRAPH_TEMPERATURE
)
from bs4 import BeautifulSoup
from PyPDF2 import PdfReader
import logging

# Configure logging
client = OpenAI(api_key=OPENAI_API_KEY)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Configure OpenAI

def modify_template(db_file, table_name, project_directory, professor_id):
    # Connect to the database
    conn = sqlite3.connect(db_file)
    conn.execute('PRAGMA foreign_keys = ON;')
    cursor = conn.cursor()

    # Fetch professor details
    cursor.execute(f"""
        SELECT "Professor", "University" FROM "{table_name}" WHERE "ID" = ?
    """, (professor_id,))
    result = cursor.fetchone()
    if not result:
        logging.error(f"No professor found with ID {professor_id}")
        conn.close()
        return

    professor_name, university = result[0], result[1]

    # Define paths
    safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
    professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
    cv_file = os.path.join(project_directory, 'Ehsan_Ghavimehr_CV.tex')

    # Ensure the professor directory exists
    if not os.path.exists(professor_dir):
        logging.error(f"Professor directory '{professor_dir}' does not exist.")
        conn.close()
        return

    # Extract and summarize text from fetched HTML pages and PDFs
    extracted_summary = extract_text_from_files(professor_dir, professor_name)

    # Save the extracted summary into a separate txt file
    notes_file = os.path.join(professor_dir, 'extracted_summary.txt')
    with open(notes_file, 'w', encoding='utf-8') as f:
        f.write(extracted_summary)
    logging.info(f"Extracted summary saved to {notes_file}")

    # Generate the personalized paragraph
    personalized_paragraph = generate_personalized_paragraph(cv_file, extracted_summary)

    # Read the email template
    template_file = os.path.join(project_directory, 'template.html')
    if not os.path.exists(template_file):
        logging.error(f"Template file '{template_file}' does not exist.")
        conn.close()
        return
    with open(template_file, 'r', encoding='utf-8') as f:
        template_content = f.read()

    # Replace placeholders
    placeholders = {
        '{{ProfessorName}}': professor_name,
        '{{University}}': university,
        '{{PersonalizedParagraph}}': personalized_paragraph,
    }
    modified_content = template_content
    for placeholder, value in placeholders.items():
        modified_content = modified_content.replace(placeholder, value)

    # Save the modified HTML in the professor's folder
    modified_email_file = os.path.join(professor_dir, 'email1.html')
    with open(modified_email_file, 'w', encoding='utf-8') as f:
        f.write(modified_content)
    logging.info(f"Modified email saved to {modified_email_file}")

    chronology_table = f"{table_name}_chronology"
    cursor.execute(f'''
        UPDATE "{chronology_table}"
        SET "html_generation_completed" = TRUE
        WHERE "ID" = ?
    ''', (professor_id,))

    # Read CV content for keywords generation
    if not os.path.exists(cv_file):
        logging.error(f"CV file '{cv_file}' does not exist.")
        conn.close()
        return

    with open(cv_file, 'r', encoding='utf-8') as f:
        cv_content = f.read()

    # Generate the prompt for CV keywords
    prompt_keywords = generate_prompt_keywords(cv_content, extracted_summary)

    # Generate the new Research Interest section
    new_research_interest = generate_research_interest(prompt_keywords)

    # Modify CV
    modify_cv(cv_file, new_research_interest, professor_dir)

    # Update the database
    cursor.execute(f'''
        UPDATE "{chronology_table}"
        SET "cv_generation_completed" = TRUE
        WHERE "ID" = ?
    ''', (professor_id,))
    conn.commit()
    conn.close()

def extract_text_from_files(professor_dir, professor_name):
    extracted_text = ''

    # Extract text from HTML files
    for filename in os.listdir(professor_dir):
        if filename.endswith('.html'):
            filepath = os.path.join(professor_dir, filename)
            with open(filepath, 'r', encoding='utf-8') as f:
                html_content = f.read()
                soup = BeautifulSoup(html_content, 'html.parser')
                text = soup.get_text(separator=' ', strip=True)
                extracted_text += text + '\n'

    # Extract text from PDF files
    for filename in os.listdir(professor_dir):
        if filename.endswith('.pdf'):
            filepath = os.path.join(professor_dir, filename)
            try:
                reader = PdfReader(filepath)
                for page in reader.pages:
                    text = page.extract_text()
                    if text:
                        extracted_text += text + '\n'
            except Exception as e:
                logging.error(f"Error extracting text from {filename}: {e}")

    # Summarize the extracted text
    summary = summarize_extracted_text(extracted_text, professor_name)
    return summary

def summarize_extracted_text(extracted_text, professor_name):
    # Use the OpenAI API to summarize the text based on professor's name
    prompt = f"""
    Professor {professor_name} has the following content extracted from their webpages and publications:

    {extracted_text}

    Please provide a concise summary (2-3 paragraphs) focusing on Professor {professor_name}'s background and key research areas.
    """

    try:
        response = client.chat.completions.create(model='GPT-3.5-turbo',
        messages=[
            {"role": "system", "content": "You are a helpful assistant that summarizes academic content."},
            {"role": "user", "content": prompt}
        ],
        max_tokens=OPENAI_SUMMARY_MAX_TOKENS,
        temperature=OPENAI_SUMMARY_TEMPERATURE)
        summary = response.choices[0].message.content.strip()
        return summary
    except Exception as e:
        logging.error(f"OpenAI API error during summarization: {e}")
        return ""

def generate_personalized_paragraph(cv_file, extracted_summary):
    # Use the OpenAI API to generate the paragraph
    with open(cv_file, 'r', encoding='utf-8') as f:
        cv_content = f.read()

    prompt = f"""
    Based on the following summarized content about the professor and my CV, write a concise paragraph explaining how my skills and background align with the professor's research interests.

    Professor's Summary:
    {extracted_summary}

    My CV:
    {cv_content}

    The paragraph should be professional and tailored to an email context.
    """

    try:
        response = client.chat.completions.create(model='gpt-4',
        messages=[
            {"role": "system", "content": "You are an assistant that writes professional emails and documents."},
            {"role": "user", "content": prompt}
        ],
        max_tokens=OPENAI_PARAGRAPH_MAX_TOKENS,
        temperature=OPENAI_PARAGRAPH_TEMPERATURE)
        generated_text = response.choices[0].message.content.strip()
        return generated_text
    except Exception as e:
        logging.error(f"OpenAI API error: {e}")
        return ""

def generate_prompt_keywords(cv_content, extracted_summary):
    prompt = f"""
    Instead of A, B, C, D, E, and F, generate 6 keywords that are shared between my background and the professor's recent research.
    My CV:
    {cv_content}

    Professor's summary:
    {extracted_summary}

    If there is low similarity between my and the professor's background, use the following keywords:
    Emotions in Decision Making, Brain Evo-Devo, Personalized Neuropsychiatry, Moral and Aesthetic Psychology, Autism, rTMS
    """

    try:
        response = client.chat.completions.create(model='gpt-4',
        messages=[
            {"role": "system", "content": "You can just answer in this structure 'A & B & C \\ D & E & F'."},
            {"role": "user", "content": prompt}
        ],
        max_tokens=35,
        temperature=0.2)
        keywords = response.choices[0].message.content.strip()
        return keywords
    except Exception as e:
        logging.error(f"OpenAI API error during keyword generation: {e}")
        return "Emotions in Decision Making & Brain Evo-Devo & Personalized Neuropsychiatry \\ Moral and Aesthetic Psychology & Autism & rTMS"

def generate_research_interest(keywords):
    try:
        prompt = f"Write a concise research interest paragraph using these keywords: {keywords}"
        response = client.chat.completions.create(model='gpt-4',
        messages=[
            {"role": "system", "content": "You are an assistant that writes professional academic content."},
            {"role": "user", "content": prompt}
        ],
        max_tokens=OPENAI_PARAGRAPH_MAX_TOKENS,
        temperature=OPENAI_PARAGRAPH_TEMPERATURE)
        research_interest = response.choices[0].message.content.strip()
        return research_interest
    except Exception as e:
        logging.error(f"OpenAI API error: {e}")
        return ""

def modify_cv(cv_file, new_research_interest, professor_dir):
    # Read the CV content
    with open(cv_file, 'r', encoding='utf-8') as f:
        cv_content = f.read()

    # Replace the old Research Interest section with the new one
    start_marker = '%BEGIN_RESEARCH_INTEREST%'
    end_marker = '%END_RESEARCH_INTEREST%'
    if start_marker in cv_content and end_marker in cv_content:
        before = cv_content.split(start_marker)[0] + start_marker + '\n'
        after = '\n' + end_marker + cv_content.split(end_marker)[1]
        cv_content = before + new_research_interest + after
    else:
        logging.error("Markers for Research Interest section not found in CV.")
        return

    # Save the modified CV in the professor's directory
    modified_cv_file = os.path.join(professor_dir, 'Ehsan_Ghavimehr_CV_modified.tex')
    with open(modified_cv_file, 'w', encoding='utf-8') as f:
        f.write(cv_content)
    logging.info(f"Modified CV saved to {modified_cv_file}")

    # Compile the modified CV using XeLaTeX
    compile_cv(modified_cv_file)

def compile_cv(tex_file):
    # Compile the TeX file using XeLaTeX
    try:
        subprocess.run(['xelatex', tex_file], cwd=os.path.dirname(tex_file), check=True)
        logging.info(f"Compiled CV saved in {os.path.dirname(tex_file)}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error compiling CV: {e}")

if __name__ == '__main__':
    # Usage: modifier.py <db_file> <table_name> <project_directory> <professor_id>
    if len(sys.argv) != 5:
        print("Usage: modifier.py <db_file> <table_name> <project_directory> <professor_id>")
        sys.exit(1)
    db_file = sys.argv[1]
    table_name = sys.argv[2]
    project_directory = sys.argv[3]
    professor_id = int(sys.argv[4])
    modify_template(db_file, table_name, project_directory, professor_id)

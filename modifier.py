# modifier.py

import os
import sys
import json
import subprocess
import openai
from openai import OpenAI
from bs4 import BeautifulSoup
from PyPDF2 import PdfReader
import logging
import sqlite3
import time  # Import time module for delays
from config import OPENAI_API_KEY

# Set OpenAI API key
client = OpenAI(api_key=OPENAI_API_KEY)

def read_simplified_cv(cv_file_path):
    with open(cv_file_path, 'r', encoding='utf-8') as f:
        cv_content = f.read()
    return cv_content

def modify_template(db_file, table_name, project_directory, professor_id):
    # Configure logging
    logger = logging.getLogger('modifier')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Remove any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler for logging
    fh = logging.FileHandler('modifier.log')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Stream handler for console output
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Connect to the database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Fetch professor details
    cursor.execute(f"""
        SELECT "Professor", "University"
        FROM "{table_name}"
        WHERE "ID" = ?
    """, (professor_id,))
    result = cursor.fetchone()
    if not result:
        logger.error(f"No professor found with ID {professor_id}")
        conn.close()
        return
    professor_name, university = result

    # Close the database connection for now
    conn.close()

    # Define paths
    safe_professor_name = ''.join(c if c.isalnum() else '_' for c in professor_name)
    professor_dir = os.path.join(project_directory, 'data', safe_professor_name)
    data_file = os.path.join(professor_dir, 'professor_data_filtered.json')
    cv_file = os.path.join(project_directory, 'Ehsan_Ghavimehr_CV.tex')  
    cv_simplified_file = os.path.join(project_directory, 'CV_simplified.txt')

    # Ensure the professor directory exists
    if not os.path.exists(professor_dir):
        logger.error(f"Professor directory '{professor_dir}' does not exist.")
        return

    # Read professor_data_filtered.json
    if not os.path.exists(data_file):
        logger.error(f"Filtered data file '{data_file}' does not exist.")
        return
    with open(data_file, 'r', encoding='utf-8') as f:
        professor_data = json.load(f)

    if not os.path.exists(cv_simplified_file):
        logger.error(f"Simplified CV file '{cv_simplified_file}' does not exist.")
        return
    simplified_cv_content = read_simplified_cv(cv_simplified_file)

    # Extract and summarize text from fetched HTML pages and PDFs
    extracted_summary = extract_text_from_files(professor_dir, logger)

    # Save the extracted summary into a separate txt file
    notes_file = os.path.join(professor_dir, 'extracted_notes.txt')
    with open(notes_file, 'w', encoding='utf-8') as f:
        f.write(extracted_summary)
    logger.info(f"Extracted summary saved to {notes_file}")

    # Generate the prompt for personalized paragraph
    prompt_paragraph = generate_prompt_paragraph(
        professor_name, professor_data, extracted_summary, simplified_cv_content
    )

    # Generate the personalized paragraph
    personalized_paragraph = generate_personalized_paragraph(prompt_paragraph, logger, model='gpt-4')

    # Read the email template
    template_file = os.path.join(project_directory, 'template.html')
    if not os.path.exists(template_file):
        logger.error(f"Template file '{template_file}' does not exist.")
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
    os.makedirs(professor_dir, exist_ok=True)
    output_file = os.path.join(professor_dir, 'email1.html')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(modified_content)
    logger.info(f"Modified email saved to {output_file}")

    # Generate the prompt for CV keywords
    prompt_keywords = generate_prompt_keywords(
        professor_name, professor_data, extracted_summary, simplified_cv_content
    )

    # Define the default keywords in case of API failure
    DEFAULT_KEYWORDS = "Emotions in Decision Making & Brain Evo-Devo & Personalized Neuropsychiatry \\\\ Moral and Aesthetic Psychology & Autism & rTMS"

    # Generate the new Research Interest section
    new_research_interest = generate_personalized_paragraph(
        prompt_keywords, logger, max_tokens=250, model='gpt-4o', default_value=DEFAULT_KEYWORDS
    )

    # Modify your CV
    modify_cv(cv_file, new_research_interest, professor_dir, logger)

    # Update the database to set html_generation and cv_generation to TRUE
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    chronology_table = f"{table_name}_chronology"
    cursor.execute(f'''
        UPDATE "{chronology_table}"
        SET "html_generation_completed" = TRUE, "cv_generation_completed" = TRUE
        WHERE "ID" = ?
    ''', (professor_id,))
    conn.commit()
    conn.close()

def extract_text_from_files(professor_dir, logger):
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

    # Extract text from PDF files using PyPDF2
    for filename in os.listdir(professor_dir):
        if filename.endswith('.pdf'):
            filepath = os.path.join(professor_dir, filename)
            try:
                reader = PdfReader(filepath)
                text = ''
                for page in reader.pages:
                    page_text = page.extract_text()
                    if page_text:
                        text += page_text + '\n'
                extracted_text += text + '\n'
            except Exception as e:
                logger.error(f"Error extracting text from {filename}: {e}")

    # Summarize the extracted text
    summary = summarize_extracted_text(extracted_text, logger)
    return summary

def summarize_extracted_text(extracted_text, logger):
    # Use the OpenAI API to summarize the text
    prompt = f"""
    Please summerize your recent research focus and current lab reseach in 3 parapheraphes. (Don't say old data about your background.)

    {extracted_text}
    """
    model = 'gpt-3.5-turbo'  # Specify the model to use
    summary = call_openai_api(prompt, max_tokens=400, model=model, role_description="You are a helpful assistant.", temperature=0.3, logger=logger)
    if summary is None:
        logger.warning("Failed to summarize extracted text after retries. Returning empty summary.")
        return ""
    return summary

def generate_prompt_paragraph(professor_name, professor_data, extracted_summary, simplified_cv_content):
    # Create a concise prompt using the simplified CV content
    prompt = f"""
    I am {professor_name}. In 3 simple sentences, explaining how your skills and background align with the my research interests. 
    Don't exagerate your skills. Don't mention the university names of research centers you have worked in.
    write nutral and humanized. Mention that you are eager to learn skills you don't yet possess. Keep it very very short.
    Don't mention neither your name nor the professor's name. Don't write neither a oppening (like Holle) nor closing (like sincerely).

    my recent articles:
    {professor_data}

    my webpages:
    {extracted_summary}

    Your CV:
    {simplified_cv_content}

    """
    return prompt

def generate_prompt_keywords(professor_name, professor_data, extracted_summary, simplified_cv_content):
    prompt = f"""
    Find 6 overlapping research topics between your background and the {professor_name}'s recent research articles (Lets call them A, B, C, D, E, and F)
    Keep keywords short

    {professor_name}'s recent articles:
    {professor_data}

    Summarized {professor_name}'s webpages:
    {extracted_summary}

    Your CV:
    {simplified_cv_content}

    seperate the 6 research topics like this format: A & B & C & D & E & F
    Don't explain anything
    If you couldn't find any research interest overlap, please use these: Emotions in Decision Making & Brain Evo-Devo & Personalized Neuropsychiatry & Moral and Aesthetic Psychology & Autism & rTMS


    """
    return prompt

def generate_personalized_paragraph(prompt, logger, max_tokens=250, model='gpt-4o', default_value=None):
    # Use the OpenAI API to generate the paragraph with retry mechanism
    response_text = call_openai_api(prompt, max_tokens=max_tokens, model=model, role_description="You are Ehsan Ghavimehr, M.D., who is applying for a Ph.D. You write simple and short.", temperature=0.7, logger=logger)
    if response_text is None:
        if default_value is not None:
            logger.warning("Failed to generate personalized paragraph after retries. Using default value.")
            return default_value
        else:
            logger.warning("Failed to generate personalized paragraph after retries. Leaving the paragraph empty.")
            return ""
    else:
        return response_text

def call_openai_api(prompt, max_tokens=200, model='gpt-4o', role_description="You are Ehsan Ghavimehr.", temperature=0.5, logger=None):
    retries = 5
    delay = 60  # Start with 60 seconds
    for attempt in range(retries):
        try:
            response = client.chat.completions.create(model=model,
            messages=[
                {"role": "system", "content": role_description},
                {"role": "user", "content": prompt}
            ],
            max_tokens=max_tokens,
            n=1,
            stop=None,
            temperature=temperature)
            generated_text = response.choices[0].message.content.strip()
            return generated_text
        except openai.RateLimitError as e:
            if logger:
                logger.warning(f"Rate limit exceeded: {e}. Retrying in {delay} seconds...")
            else:
                print(f"Rate limit exceeded: {e}. Retrying in {delay} seconds...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
        except Exception as e:
            if logger:
                logger.error(f"OpenAI API error: {e}")
            else:
                print(f"OpenAI API error: {e}")
            return None
    if logger:
        logger.error("Failed to get a response from OpenAI API after multiple retries.")
    else:
        print("Failed to get a response from OpenAI API after multiple retries.")
    return None

def modify_cv(cv_file, new_research_interest, professor_dir, logger):
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
        logger.warning("Markers for Research Interest section not found in CV.")
        return

    # Save the modified CV in the professor's directory
    modified_cv_file = os.path.join(professor_dir, 'Ehsan_Ghavimehr_CV.tex')
    with open(modified_cv_file, 'w', encoding='utf-8') as f:
        f.write(cv_content)
    logger.info(f"Modified CV saved to {modified_cv_file}")

    # Compile the modified CV using XeLaTeX
    compile_cv(modified_cv_file, logger)

def compile_cv(tex_file, logger):
    # Compile the TeX file using XeLaTeX
    try:
        subprocess.run(['xelatex', '-interaction=nonstopmode', tex_file], cwd=os.path.dirname(tex_file), check=True)
        logger.info(f"Compiled CV saved in {os.path.dirname(tex_file)}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error compiling CV: {e}")

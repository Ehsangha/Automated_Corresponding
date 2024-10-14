# config.py

# API Keys
SERPAPI_API_KEY = 'your_key'
LINKEDIN_EMAIL = 'you'
LINKEDIN_PASSWORD = 'pass'
OPENAI_API_KEY = 'your_key'
ORCID_CLIENT_ID = 'your_key'
ORCID_CLIENT_SECRET = 'your_key'
ENTREZ_EMAIL = 'your_email'


# Database and Project Configurations
DB_FILE = 'main.db'                     # Default database file
TABLE_NAME = 'neuro'      # Default table name
SEARCH_DEPTH = 0                         # Default search depth
PROJECT_DIRECTORY = 'path/to/automated_corresponding'  # Default project directory

# Search Style (if applicable)
SEARCH_STYLE = 1           # or 'depth-first'

# Reminder Intervals
REMINDER_INTERVAL_1 = 7    # days until first reminder
REMINDER_INTERVAL_2 = 14   # days until second reminder
REMINDER_INTERVAL_3 = 21   # days until third reminder

# Email Settings
EMAIL_ACCOUNTS = [
    {
        "from_email": 'your_email',
        "username": 'your_email',
        "password": 'your_pass',
        "smtp_host": 'smtp.google.com',
        "smtp_port": 465,
        "imap_host": 'imap.google.com', 
        "imap_port": 993,
        "ssl": True,
    },
    {
        "from_email": 'your_email',
        "username": 'your_email',
        "password": 'your_pass',
        "smtp_host": 'smtp.google.com',
        "smtp_port": 465,
        "imap_host": 'imap.google.com', 
        "imap_port": 993,
        "ssl": True,
    },
    {
        "from_email": 'your_email',
        "username": 'your_email',
        "password": 'your_pass',
        "smtp_host": 'smtp.provider.com',
        "smtp_port": 465,
        "imap_host": 'imap.provider.com', 
        "imap_port": 993,
        "ssl": True,
    },
]

# Test Run Configuration
TEST_RUN = True  # If True, emails will be sent to TEST_EMAIL instead of the professor's real emails
TEST_EMAIL = 'your_email'

SENDING_METHOD = 'smtp'  


# Temperature and max tokens for email generation
EMAIL_PERSONALIZED_PARAGRAPH_TEMPERATURE = 0.3
EMAIL_PERSONALIZED_PARAGRAPH_MAX_TOKENS = 250

# Temperature and max tokens for CV modification
CV_RESEARCH_INTEREST_TEMPERATURE = 0.3
CV_RESEARCH_INTEREST_MAX_TOKENS = 100

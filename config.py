# config.py

# API Keys
SERPAPI_API_KEY = '976d474c8af3e91139e205333918849e95010795cb5535631f1ec0edd6de1f48'
LINKEDIN_EMAIL = 'ghavimehr'
LINKEDIN_PASSWORD = ';hlfd,jv'
OPENAI_API_KEY = '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'
ORCID_CLIENT_ID = 'APP-PPC15PE7WPU5RMYL'
ORCID_CLIENT_SECRET = 'e3f50710-712c-43ac-8a70-94cfe55c4434'
ENTREZ_EMAIL = 'ehsanghavimehr@gmail.com'


# Database and Project Configurations
DB_FILE = 'main.db'                     # Default database file
TABLE_NAME = 'neuro'      # Default table name
SEARCH_DEPTH = 0                         # Default search depth
PROJECT_DIRECTORY = '/home/ehsan/usr/automated_corresponding'  # Default project directory

# Search Style (if applicable)
SEARCH_STYLE = 1           # or 'depth-first'

# Reminder Intervals
REMINDER_INTERVAL_1 = 7    # days until first reminder
REMINDER_INTERVAL_2 = 14   # days until second reminder
REMINDER_INTERVAL_3 = 21   # days until third reminder

# Email Settings
EMAIL_ACCOUNTS = [
    {
        "from_email": 'eh.ghavimehr@gmail.com',
        "username": 'eh.ghavimehr@gmail.com',
        "password": 'dbugdhtzsrvqinpp',
        "smtp_host": 'smtp.google.com',
        "smtp_port": 465,
        "imap_host": 'imap.google.com', 
        "imap_port": 993,
        "ssl": True,
    },
    {
        "from_email": 'E.Ghavimehr@gmail.com',
        "username": 'E.Ghavimehr@gmail.com',
        "password": 'kmaloqwfjdxtnjxi',
        "smtp_host": 'smtp.google.com',
        "smtp_port": 465,
        "imap_host": 'imap.google.com', 
        "imap_port": 993,
        "ssl": True,
    },
    {
        "from_email": 'ehsan@ghavimehr.com',
        "username": 'ehsan@ghavimehr.com',
        "password": 'Now;hlfd,jv1',
        "smtp_host": 'smtp.hostinger.com',
        "smtp_port": 465,
        "imap_host": 'imap.hostinger.com', 
        "imap_port": 993,
        "ssl": True,
    },
]

# Test Run Configuration
TEST_RUN = False  # If True, emails will be sent to TEST_EMAIL instead of the professor's real emails
TEST_EMAIL = 'ehsanghavimehr@gmail.com'

SENDING_METHOD = 'smtp'  


# Temperature and max tokens for email generation
EMAIL_PERSONALIZED_PARAGRAPH_TEMPERATURE = 0.3
EMAIL_PERSONALIZED_PARAGRAPH_MAX_TOKENS = 250

# Temperature and max tokens for CV modification
CV_RESEARCH_INTEREST_TEMPERATURE = 0.3
CV_RESEARCH_INTEREST_MAX_TOKENS = 100
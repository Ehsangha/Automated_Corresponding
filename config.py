#!/home/ehsan/anaconda3/bin/python3


# API Keys
SERPAPI_API_KEY = '0fc3a65cc2b26bb9645cba8a589df5ac50f926cfc9a1719a7dc3478b65c88205'
LINKEDIN_EMAIL = ''
LINKEDIN_PASSWORD = ''

OPENAI_API_KEY = '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'

ORCID_CLIENT_ID = 'APP-PPC15PE7WPU5RMYL'
ORCID_CLIENT_SECRET = 'e3f50710-712c-43ac-8a70-94cfe55c4434'

# Other configurations
ENTREZ_EMAIL = 'ehsanghavimehr@gmail.com'

# Database and Project Configurations
DB_FILE = 'main.db'                     # Default database file
TABLE_NAME = 'cognitive_psychology'      # Default table name
SEARCH_DEPTH = 1                         # Default search depth
PROJECT_DIRECTORY = '/home/ehsan/usr/mailling_list'  # Default project directory

# Search Style (if applicable)
SEARCH_STYLE = 'breadth-first'           # or 'depth-first'

# Reminder Intervals
REMINDER_INTERVAL_1 = 7    # days until first reminder
REMINDER_INTERVAL_2 = 14   # days until second reminder
REMINDER_INTERVAL_3 = 21   # days until third reminder

# Email Settings
EMAIL_ACCOUNTS = [
    {
        'from_email': 'ehsan@ghavimehr.com',
        'username': 'ehsan@ghavimehr.com',
        'password': '',
        'host': 'smtp.hostinger.com',  # Hostinger's SMTP server
        'port': 465,  # Typically 465 for SSL
        'use_ssl': True,
        'imap_host': 'imap.hostinger.com',  # Hostinger's IMAP server
        'imap_port': 993,  # Typically 993 for SSL
    },
    # You can add more email accounts here
]

SENDING_METHOD = 1

# Test Run Configuration
TEST_RUN = True  # If True, emails will be sent to TEST_EMAIL instead of the professor's real emails
TEST_EMAIL = 'ehsanghavimehr@gmail.com'

# OpenAI Configuration
# OpenAI Configuration for Summarization
OPENAI_SUMMARY_MAX_TOKENS = 1000
OPENAI_SUMMARY_TEMPERATURE = 0.3

# OpenAI Configuration for Paragraph Generation
OPENAI_PARAGRAPH_MAX_TOKENS = 125
OPENAI_PARAGRAPH_TEMPERATURE = 0.5

# OpenAI Configuration for Keywords Generation
OPENAI_KEYWORDS_MAX_TOKENS = 35
OPENAI_KEYWORDS_TEMPERATURE = 0.2

# Answer State (assuming initial state is 0)
ANSWER_STATE = 0  # 0: 'No response yet'

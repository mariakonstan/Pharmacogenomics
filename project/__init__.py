import logging
from .settings import LOG_FILE


# Create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)  # Set lowest level to capture everything

# Create formatter
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Create file handler for logging to file
file_handler = logging.FileHandler(LOG_FILE)
file_handler.setLevel(logging.DEBUG)  # Log everything to file
file_handler.setFormatter(formatter)

# Create console handler for logging to console (warnings and errors only)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)  # Only show WARNING and above in console
console_handler.setFormatter(formatter)

# Add handlers to the logger (avoid duplicate handlers)
if not logger.handlers:
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logger.info("Logger initialized.")

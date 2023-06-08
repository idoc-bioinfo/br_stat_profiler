import logging
import time
import sys
from humanize import precisedelta

# Create a custom log formatter
class ElapsedTimeFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.previous_time = None
        self.start_time = time.perf_counter()
        
    def format(self, record):
        current_time = time.perf_counter()
        elapsed_msg = "(  Init   )"  # first msg
        secs_from_start = int(current_time - self.start_time)
        mins_from_start = int(secs_from_start /60)
        secs_reminder = secs_from_start - (mins_from_start * 60)
        timer_msg = f"{mins_from_start}m:{secs_reminder}s"
        
        if self.previous_time is not None:
            elapsed_time = current_time - self.previous_time
            elapsed_msg = f"(+{elapsed_time:.4f} s)"
            
        self.previous_time = current_time
        # Format the log message with the elapsed time
        formatted_message = super().format(record)
        formatted_message = f'{timer_msg} {elapsed_msg} - {formatted_message}'
    
        return formatted_message

logging.getLogger().handlers = []
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.disable(logging.CRITICAL)


def initialize_logger(file_descriptor=sys.stderr, level=logging.INFO):
    logging.disable(level-1)
    # Set the custom formatter for the file handler
    file_formatter = ElapsedTimeFormatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.StreamHandler(file_descriptor)
    file_handler.setFormatter(file_formatter)
    # Add the file handler to the logger
    logger.addHandler(file_handler)
    
def get_logger():
    return logging.getLogger(__name__)

  
if __name__ == "__main__":

    # Create a file handler with an open file descriptor
    # file_descriptor = open('log.txt', 'a')  # Replace with your desired file name and mode

    initialize_logger(level=logging.DEBUG)

    # initialize_logger(file_descriptor)
    # file_handler = logging.StreamHandler(file_descriptor)
    # file_handler.setFormatter(file_formatter)

# Log messages with previous time
    logger.debug('This is the first log message')
    time.sleep(1)  # Simulating delay

    logger.info('This is the second log message')
    time.sleep(1)  # Simulating delay
    # logging.disable(logging.CRITICAL)

    logger.debug('This is the third log message')
    time.sleep(1)  # Simulating delay

    logger.warning('This is the forth log message')

import logging
import sys


class RedirectStdoutToLogger:
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self._original_stdout = sys.stdout

    def write(self, message):
        # Log the message only if it contains non-whitespace characters
        if message.strip():  # Check if the message is not empty
            self.logger.log(
                self.log_level, message.strip()
            )  # Log without removing leading/trailing spaces

    def flush(self):
        pass  # No need to flush anything

    def __enter__(self):
        sys.stdout = self  # Redirect stdout to this object
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self._original_stdout  # Reset stdout when done

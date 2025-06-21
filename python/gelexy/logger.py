import logging
import sys

level = logging.INFO


class LowercaseLevelnameFormatter(logging.Formatter):
    def format(self, record):
        GREEN = "\033[92m"
        RESET = "\033[0m"

        record.levelname = f"{GREEN}{record.levelname.lower()}{RESET}"
        return super().format(record)


def setup_logger(name: str):
    logger = logging.getLogger(name)

    if not logger.hasHandlers():
        handler = logging.StreamHandler(sys.stdout)
        if level == logging.INFO:
            handler.setFormatter(LowercaseLevelnameFormatter("%(message)s"))
        else:
            handler.setFormatter(
                LowercaseLevelnameFormatter(
                    "%(asctime)s %(levelname)s %(message)s"
                )
            )
        logger.addHandler(handler)
        logger.setLevel(level)

    return logger

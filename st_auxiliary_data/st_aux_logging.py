"""Provides special formatting for log messages
"""

import logging


class LoggingFilter(logging.Filter):
    """Forces 'ESPA-AUX' to be provided in the 'subsystem' tag of the log
       format string
    """

    def filter(self, record):
        """Provide the string for the 'subsystem' tag"""

        record.subsystem = 'ESPA-AUX'

        return True


class ExceptionFormatter(logging.Formatter):
    """Modifies how exceptions are formatted
    """

    def formatException(self, exc_info):
        """Specifies how to format the exception text"""

        result = super(ExceptionFormatter, self).formatException(exc_info)

        return repr(result)

    def format(self, record):
        """Specifies how to format the message text if it is an exception"""

        s = super(ExceptionFormatter, self).format(record)
        if record.exc_text:
            s = s.replace('\n', ' ')
            s = s.replace('\\n', ' ')

        return s

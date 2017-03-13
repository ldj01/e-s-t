"""Provides a simplified Http Session class
"""

import logging
import requests
from contextlib import closing


logger = logging.getLogger(__name__)


class HttpSession(object):
    """Manages an http(s) session.
    """

    def __init__(self, max_retries=3, block_size=None, timeout=300.0):
        """Initialization of the class

        Args:
            max_retries <int>: The number of max retries
            block_size <int>: How many bytes to read at a time
            timeout <float>: How long to wait before trying again
        """
        super(HttpSession, self).__init__()

        self.session = requests.Session()

        self.timeout = timeout
        self.max_retries = max_retries
        self.status_code = requests.codes['ok']

        # Determine if we are streaming or not based on block_size
        self.block_size = block_size
        self.stream = False
        if self.block_size is not None:
            self.stream = True

        adapter = requests.adapters.HTTPAdapter(max_retries=max_retries)
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)

    def login(self, login_url, login_data):
        """Provides for establishing a logged in session.

        Args:
            login_url <str>: The URL to use for logging in to the system
            login_data <dict>: Dictionary of the login data fields
        """

        # Login to the site
        self.session.post(url=login_url, data=login_data)

    def _get_file(self, download_url, destination_file, headers=None):
        """Download the whole source file into memory before dumping to
           the local file

        Args:
            download_url <str>: The URL to download
            destination_file <str>: The file to create on the localhost
            headers <str>: Any special headers to pass to the session
        """

        with closing(self.session.get(url=download_url,
                                      timeout=self.timeout,
                                      headers=headers)) as req:

            self.status_code = req.status_code

            if not req.ok:
                logger.error('HTTP - Transfer of [{0}] - FAILED'
                             .format(download_url))
                # The raise_for_status generates an exception to be caught
                req.raise_for_status()

            # Write the downloaded data to the destination file
            with open(destination_file, 'wb') as local_fd:
                local_fd.write(req.content)

    def _stream_file(self, download_url, destination_file, headers=None):
        """Download by streaming 'block_size' of data at a time

        Args:
            download_url <str>: The URL to download
            destination_file <str>: The file to create on the localhost
            headers <str>: Any special headers to pass to the session
        """

        retrieved_bytes = 0
        with closing(self.session.get(url=download_url,
                                      timeout=self.timeout,
                                      stream=True,
                                      headers=headers)) as req:
            self.status_code = req.status_code

            if not req.ok:
                logger.error('HTTP - Transfer of [{0}] - FAILED'
                             .format(download_url))
                # The raise_for_status generates an exception to be caught
                req.raise_for_status()

            file_size = int(req.headers['content-length'])

            # Set block size based on streaming
            if self.stream:
                block_size = self.block_size
            else:
                block_size = file_size

            # Write the downloaded data to the destination file
            with open(destination_file, 'wb') as local_fd:
                for data_chunk in req.iter_content(block_size):
                    local_fd.write(data_chunk)
                    retrieved_bytes += len(data_chunk)

            if retrieved_bytes != file_size:
                raise Exception('Transfer Failed - HTTP -'
                                ' Retrieved {0} out of {1} bytes'
                                .format(retrieved_bytes, file_size))

    def http_transfer_file(self, download_url, destination_file):
        """Use http to transfer a file from a source location to a
           destination file on the localhost

        Args:
            download_url <str>: The URL to download
            destination_file <str>: The file to create on the localhost

        Returns:
            status_code - One of the following
                        - 200, requests.codes['ok']
                        - 404, requests.codes['not_found']:
                        - 503, requests.codes['service_unavailable']:

        Notes:
            If a 503 is returned, the logged exception should be reviewed
            to determine the real cause of the error.
        """

        logger.info(download_url)

        retry_attempt = 0
        done = False
        while not done:
            self.status_code = requests.codes['ok']
            try:

                self._stream_file(download_url, destination_file)

                logger.info("Transfer Complete - HTTP")
                done = True

            except IOError:
                logger.exception('HTTP - Transfer Issue')

                if self.status_code not in (requests.codes['not_found'],
                                            requests.codes['forbidden']):

                    if retry_attempt > self.max_retries:
                        logger.info('HTTP - Transfer Failed'
                                    ' - exceeded retry limit')
                        done = True
                    else:
                        retry_attempt += 1
                        sleep(int(1.5 * retry_attempt))
                else:
                    # Not Found - So break the looping because we are done
                    done = True

        return self.status_code

    def get_lines_from_url(self, download_url):
        """Retrieve lines from a url

        Args:
            download_url <str>: The URL to download

        Returns:
            <list>: A list containing the lines from the URL
        """

        data = list()
        self.status_code = requests.codes['ok']

        with closing(self.session.get(url=download_url,
                                      timeout=self.timeout,
                                      stream=self.stream)) as req:
            self.status_code = req.status_code

            if not req.ok:
                logger.error('HTTP - Transfer of [{0}] - FAILED'
                             .format(download_url))
                # The raise_for_status generates an exception to be caught
                req.raise_for_status()

            for line in req.iter_lines():
                data.append(line)

        return data

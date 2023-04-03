import logging

pyvbrc_log = logging.getLogger("pyVBRc")
pyvbrc_log.setLevel(logging.INFO)

_formatter = logging.Formatter("%(name)s : [%(levelname)s ] %(asctime)s:  %(message)s")

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(_formatter)
pyvbrc_log.addHandler(stream_handler)

import logging
import os

CATHPY_DEBUG=os.environ['CATHPY_DEBUG'] if 'CATHPY_DEBUG' in os.environ else False 

#logging.basicConfig(level=logging.INFO)

debug = logging.StreamHandler()
debug.setLevel(logging.DEBUG)
debug_format = logging.Formatter('%(asctime)s %(name)15s:%(lineno)-5s %(levelname)-8s | %(message)s')
debug.setFormatter(debug_format)

info = logging.StreamHandler()
info.setLevel(logging.INFO)
info_format = logging.Formatter('%(asctime)s %(levelname)-8s | %(message)s',
    datefmt='%Y-%m-%d %H:%M')
info.setFormatter(info_format)

# add the handler to the root logger
root = logging.getLogger('')
if CATHPY_DEBUG:
    root.setLevel(logging.DEBUG)
    root.addHandler(debug)
else:
    root.setLevel(logging.INFO)
    root.addHandler(info)

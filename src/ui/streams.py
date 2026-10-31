"""Keep terminal encoding; safely replace unsupported glyphs in legacy pipes."""
import sys

def configure_streams():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")

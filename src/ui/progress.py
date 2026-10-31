"""Progress advances only when work finishes; no synthetic Vina percentages."""
from functools import wraps
from rich.progress import Progress, TextColumn, BarColumn, MofNCompleteColumn, TimeElapsedColumn
from src.ui.console import console

class WorkProgress:
    def __init__(self, description, total=None):
        self.progress = Progress(TextColumn("{task.description}"), BarColumn(),
            MofNCompleteColumn(), TimeElapsedColumn(), console=console)
        self.task = self.progress.add_task(description, total=total)

    def __enter__(self):
        self.progress.start()
        return self

    def update(self, **values):
        self.progress.update(self.task, **values)

    def advance(self):
        self.progress.advance(self.task)

    def __exit__(self, *error):
        self.progress.stop()

def evolution_progress(method):
    @wraps(method)
    def run(self, *args, **kwargs):
        with WorkProgress("Evolution: initialization") as progress:
            self._progress = progress
            try:
                return method(self, *args, **kwargs)
            finally:
                self._progress = None
    return run

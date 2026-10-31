"""console - FragEvoLab module."""


from rich.console import Console
from rich.text import Text

from src.ui.streams import configure_streams
configure_streams()

# Global Rich console
console = Console()

def print_header(title: str, subtitle: str = None, color: str = "cyan"):
    """
    Print a formatted header.
    
    Args:
        title: Title
        subtitle: Subtitle
        color: Header color
    """
    console.print(f"\n[bold {color}]{title}[/bold {color}]")
    if subtitle:
        console.print(subtitle, style="dim", markup=False)

def print_success(message: str):
    """
    Print a success message.
    
    Args:
        message: Message
    """
    console.print(f"[green]✓[/green] {message}")

def print_error(message: str):
    """
    Print an error message.
    
    Args:
        message: Message
    """
    console.print(f"[red]✗[/red] {message}")

def print_warning(message: str):
    """
    Print a warning.
    
    Args:
        message: Message
    """
    console.print(f"[yellow]⚠[/yellow] {message}")

def print_info(message: str):
    """
    Print an informational message.
    
    Args:
        message: Message
    """
    console.print(f"[blue]ℹ[/blue] {message}")


def print_app_header():
    """Print the application banner at startup."""
    from src.version import APP_NAME, __version__
    glyphs = {
        "F": ["#####", "#    ", "#### ", "#    ", "#    "],
        "R": ["#### ", "#   #", "#### ", "#  # ", "#   #"],
        "A": [" ### ", "#   #", "#####", "#   #", "#   #"],
        "G": [" ####", "#    ", "# ###", "#   #", " ####"],
        "E": ["#####", "#    ", "#### ", "#    ", "#####"],
        "V": ["#   #", "#   #", "#   #", " # # ", "  #  "],
        "O": [" ### ", "#   #", "#   #", "#   #", " ### "],
        "L": ["#    ", "#    ", "#    ", "#    ", "#####"],
        "B": ["#### ", "#   #", "#### ", "#   #", "#### "],
    }
    # ASCII only: no dependence on special molecular glyphs or terminal fonts.
    emblem = [r"     __            __",
              r" R--/  \    ~     /  \ ",
              r"    |()| -- ~ -- |()|",
              r"    \__/    ~     \__/--R1"]
    width = min(console.width, 59)
    if console.width >= 32:
        padding = " " * max(0, (width - max(map(len, emblem))) // 2)
        for line in emblem:
            console.print(padding + line, style="cyan", markup=False, highlight=False)
    if console.width >= 60:
        console.print()
        for row in range(5):
            text = Text()
            for index, letter in enumerate("FRAGEVOLAB"):
                color = "cyan" if index < 4 else "yellow" if index < 7 else "magenta"
                text.append(glyphs[letter][row], style=color)
                if index < 9: text.append(" ")
            console.print(text)
        console.print()
    console.print(f"{APP_NAME} {__version__}", style="bold cyan", markup=False)
    console.print("Evolutionary generation and docking", style="dim")

def print_section_divider(title: str = None):
    """
    Print a section separator.
    
    Args:
        title: Section title
    """
    if title:
        console.print(title, style="bold cyan", markup=False)


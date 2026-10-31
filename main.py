"""FragEvoLab launcher; --check runs diagnostics without starting calculations."""
import argparse
from pathlib import Path
from src.version import APP_NAME, __version__


def main():
    from src.ui.streams import configure_streams
    configure_streams()
    parser = argparse.ArgumentParser(description=f'{APP_NAME} {__version__}')
    parser.add_argument('--version', action='version', version=f'{APP_NAME} {__version__}')
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--check', action='store_true', help='Check installation and evolution inputs without running calculations')
    parser.add_argument('--config', type=Path, help='Settings JSON for --check')
    parser.add_argument('--seeds', type=Path, help='Seed molecule CSV for --check')
    parser.add_argument('--resume', type=Path, help='Checkpoint for --check')
    args = parser.parse_args()
    from src.core.preflight import check_project, check_dependencies, display
    if args.check:
        report = check_project(args.project, config_path=args.config, seeds=args.seeds, resume=args.resume)
        display(report)
        return 0 if report.ready else 1
    report = check_dependencies()
    if not report.ready:
        display(report)
        return 1
    from src.ui.application import FragEvoLabApp
    app = None
    try:
        app = FragEvoLabApp(args.project)
        app.run()
        return 0
    except KeyboardInterrupt:
        print('Program stopped by the user.')
        return 130
    except Exception as error:
        print(f'Startup error: {error}. Run python main.py --check for diagnostics')
        return 1
    finally:
        if app is not None:
            app.close()


if __name__ == '__main__':
    raise SystemExit(main())

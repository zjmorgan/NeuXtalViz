import argparse

from NeuXtalViz.qt.gui import gui
from NeuXtalViz.trame.gui import trame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "frontend",
        choices=["qt", "trame"],
        help="The frontend in which to display the application.",
    )

    args = parser.parse_args()

    match args.frontend:
        case "qt":
            gui()
        case "trame":
            trame()


if __name__ == "__main__":
    main()

import sys
import os
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QFileDialog,
    QTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
    QListWidget,
    QLabel,
    QHBoxLayout,
)
import subprocess


import qdarktheme

qdarktheme.enable_hi_dpi()


class CommandBrowser(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Command Browser")
        self.setGeometry(100, 100, 800, 600)

        self.directory_path = None
        self.current_file_path = None

        self.init_ui()

    def init_ui(self):
        central_widget = QWidget()
        layout = QHBoxLayout()

        self.file_list = QListWidget()
        self.file_list.itemClicked.connect(self.load_file)
        layout.addWidget(self.file_list, 1)

        editor_layout = QVBoxLayout()
        self.file_label = QLabel("No file selected")
        self.text_editor = QTextEdit()
        editor_layout.addWidget(self.file_label)
        editor_layout.addWidget(self.text_editor, 4)

        button_layout = QHBoxLayout()
        save_button = QPushButton("Save File")
        save_button.clicked.connect(self.save_file)
        self.run_button = QPushButton("Run Command")
        self.run_button.clicked.connect(self.run_command)
        button_layout.addWidget(save_button)
        button_layout.addWidget(self.run_button)
        editor_layout.addLayout(button_layout)

        layout.addLayout(editor_layout, 3)

        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        menu = self.menuBar()
        file_menu = menu.addMenu("File")
        open_dir_action = file_menu.addAction("Open Directory")
        open_dir_action.triggered.connect(self.open_directory)
        command_menu = menu.addMenu("Command")
        xprep_action = command_menu.addAction("xprep")
        xprep_action.triggered.connect(self.switch_xprep)
        shelxl_action = command_menu.addAction("shelxl")
        shelxl_action.triggered.connect(self.switch_shelxl)
        shelxt_action = command_menu.addAction("shelxt")
        shelxt_action.triggered.connect(self.switch_shelxt)
        discus_action = command_menu.addAction("discus")
        discus_action.triggered.connect(self.switch_discus)

        self.switch_xprep()

    def switch_command(self, command):
        self.command = command
        self.run_button.setText(self.command)

    def switch_xprep(self):
        self.switch_command("xprep")
        self.terminal = "{} {}"

    def switch_shelxl(self):
        self.switch_command("shelxl")
        self.terminal = "{} {}"

    def switch_shelxt(self):
        self.switch_command("shelxt")
        self.terminal = "{} {}"

    def switch_discus(self):
        self.switch_command("discus_suite")
        self.terminal = "{}"

    def open_directory(self):
        self.directory_path = QFileDialog.getExistingDirectory(
            self, "Select Directory"
        )
        if self.directory_path:
            self.load_files_in_directory()

    def load_files_in_directory(self):
        self.file_list.clear()
        for file_name in os.listdir(self.directory_path):
            self.file_list.addItem(file_name)

    def load_file(self, item):
        self.current_file_path = os.path.join(self.directory_path, item.text())
        if os.path.isfile(self.current_file_path):
            with open(self.current_file_path, "r") as file:
                self.text_editor.setText(file.read())
            self.file_label.setText(f"Editing: {item.text()}")

    def save_file(self):
        if self.current_file_path:
            with open(self.current_file_path, "w") as file:
                file.write(self.text_editor.toPlainText())
            self.statusBar().showMessage(
                f"File {os.path.basename(self.current_file_path)} saved successfully!",
                3000,
            )

    def run_command(self):
        if self.current_file_path:
            command = self.terminal.format(
                self.command, self.current_file_path
            )

            terminal_commands = [
                f"gnome-terminal -- bash -c '{command}; exec bash'",
                f"konsole -e bash -c '{command}; exec bash'",
                f"xterm -hold -e '{command}'",
                f"lxterminal -e '{command}'",
                f"xfce4-terminal -e '{command}'",
            ]

            for terminal_cmd in terminal_commands:
                try:
                    process = subprocess.Popen(terminal_cmd, shell=True)
                    if process is not None:
                        return
                except FileNotFoundError:
                    continue

            self.statusBar().showMessage("No suitable terminal found!", 5000)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    qdarktheme.setup_theme("light")
    editor = CommandBrowser()
    editor.show()
    sys.exit(app.exec_())

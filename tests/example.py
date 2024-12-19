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
from PyQt5.QtCore import Qt
from subprocess import Popen, PIPE


class FileEditorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("File Editor")
        self.setGeometry(100, 100, 800, 600)

        self.directory_path = None
        self.current_file_path = None

        self.init_ui()

    def init_ui(self):
        # Layout setup
        central_widget = QWidget()
        layout = QHBoxLayout()

        # File list
        self.file_list = QListWidget()
        self.file_list.itemClicked.connect(self.load_file)
        layout.addWidget(self.file_list, 1)

        # File editor
        editor_layout = QVBoxLayout()
        self.file_label = QLabel("No file selected")
        self.text_editor = QTextEdit()
        editor_layout.addWidget(self.file_label)
        editor_layout.addWidget(self.text_editor, 4)

        # Buttons
        button_layout = QHBoxLayout()
        save_button = QPushButton("Save File")
        save_button.clicked.connect(self.save_file)
        run_button = QPushButton("Run Command")
        run_button.clicked.connect(self.run_command)
        button_layout.addWidget(save_button)
        button_layout.addWidget(run_button)
        editor_layout.addLayout(button_layout)

        layout.addLayout(editor_layout, 3)

        # Set central widget
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        # Menu
        menu = self.menuBar()
        file_menu = menu.addMenu("File")
        open_dir_action = file_menu.addAction("Open Directory")
        open_dir_action.triggered.connect(self.open_directory)

    def open_directory(self):
        self.directory_path = QFileDialog.getExistingDirectory(
            self, "Select Directory"
        )
        if self.directory_path:
            self.load_files_in_directory()

    def load_files_in_directory(self):
        self.file_list.clear()
        for file_name in os.listdir(self.directory_path):
            # if file_name.endswith(".txt"):  # Filter for text files
            self.file_list.addItem(file_name)

    def load_file(self, item):
        self.current_file_path = os.path.join(self.directory_path, item.text())
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
            command = f"cat {self.current_file_path}"  # Example: Change this to your actual command
            process = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()

            # Display command output
            if process.returncode == 0:
                self.text_editor.setText(stdout.decode())
                self.statusBar().showMessage(
                    "Command executed successfully!", 3000
                )
            else:
                self.text_editor.setText(stderr.decode())
                self.statusBar().showMessage("Command execution failed.", 3000)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    editor = FileEditorApp()
    editor.show()
    sys.exit(app.exec_())

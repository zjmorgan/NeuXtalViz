import sys
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QTextEdit,
    QLineEdit,
    QPushButton,
)
from PyQt5.QtCore import QProcess


class BashShellWidget(QWidget):
    def __init__(self):
        super().__init__()

        # Set up the UI
        self.setWindowTitle("Bash Shell in PyQt5")
        self.layout = QVBoxLayout()

        self.output_display = QTextEdit()
        self.output_display.setReadOnly(True)
        self.layout.addWidget(self.output_display)

        self.command_input = QLineEdit()
        self.command_input.setPlaceholderText("Enter a Bash command...")
        self.layout.addWidget(self.command_input)

        self.run_button = QPushButton("Run Command")
        self.layout.addWidget(self.run_button)

        self.setLayout(self.layout)

        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.on_ready_read_output)
        self.process.readyReadStandardError.connect(self.on_ready_read_error)

        self.run_button.clicked.connect(self.run_command)
        self.command_input.returnPressed.connect(self.run_command)

    def run_command(self):
        command = self.command_input.text().strip()
        if command:
            self.output_display.append(
                f"$ {command}"
            )  # Display the command in the output window
            self.process.start(
                "bash", ["-c", command]
            )  # Execute the command in Bash
            self.command_input.clear()

    def on_ready_read_output(self):
        output = self.process.readAllStandardOutput().data().decode()
        self.output_display.append(output)

    def on_ready_read_error(self):
        error = self.process.readAllStandardError().data().decode()
        self.output_display.append(f"Error: {error}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    shell = BashShellWidget()
    shell.show()
    sys.exit(app.exec_())

stylesheet = """
/* General settings for the whole application */
QWidget {
    background-color: #f0f0f0;
    color: #333;
    font-family: Arial, Helvetica, sans-serif;
    font-size: 14px;
}

/* Style for QPushButton */
QPushButton {
    background-color: #4CAF50;
    border: none;
    color: white;
    padding: 10px 20px;
    text-align: center;
    text-decoration: none;
    font-size: 14px;
    margin: 4px 2px;
    border-radius: 12px;
}

QPushButton:hover {
    background-color: #45a049;
}

QPushButton:pressed {
    background-color: #388e3c;
}

/* Style for QLineEdit */
QLineEdit {
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 8px;
    background-color: #fff;
    color: #333;
}

QLineEdit:focus {
    border: 1px solid #4CAF50;
}

/* Style for QLabel */
QLabel {
    font-size: 14px;
    color: #333;
}

/* Style for QMenuBar */
QMenuBar {
    background-color: #333;
    color: white;
    padding: 5px;
}

QMenuBar::item {
    background-color: #333;
    color: white;
}

QMenuBar::item:selected {
    background-color: #4CAF50;
}

/* Style for QStatusBar */
QStatusBar {
    background-color: #333;
    color: white;
    padding: 5px;
}

/* Style for QComboBox */
QComboBox {
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 8px;
    background-color: #fff;
    color: #333;
}

QComboBox:focus {
    border: 1px solid #4CAF50;
}

QComboBox QAbstractItemView {
    border: 1px solid #4CAF50;
    selection-background-color: #4CAF50;
    selection-color: white;
}

/* Style for QTabWidget */
QTabWidget::pane {
    border: 1px solid #ccc;
    padding: 5px;
}

QTabWidget::tab-bar {
    left: 5px;
}

QTabBar::tab {
    background: #ccc;
    border: 1px solid #ccc;
    padding: 8px;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    color: #333;
}

QTabBar::tab:selected, QTabBar::tab:hover {
    background: #4CAF50;
    color: white;
}

/* Scrollbar */
QScrollBar:vertical {
    border: 1px solid #ccc;
    background: #fff;
    width: 12px;
    margin: 0px 0px 0px 0px;
}

QScrollBar::handle:vertical {
    background: #4CAF50;
    min-height: 20px;
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    background: #ccc;
    height: 0px;
}

QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
    background: none;
}

QScrollBar:horizontal {
    border: 1px solid #ccc;
    background: #fff;
    height: 12px;
    margin: 0px 0px 0px 0px;
}

QScrollBar::handle:horizontal {
    background: #4CAF50;
    min-width: 20px;
}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    background: #ccc;
    width: 0px;
}

QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {
    background: none;
}
"""
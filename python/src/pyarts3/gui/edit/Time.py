"""Editor for Time type."""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QDialogButtonBox, QDateTimeEdit
from PyQt5.QtCore import QDateTime
import datetime


def edit(value, parent=None):
    """
    Edit a Time object using a datetime picker.
    
    The Time object has a .time attribute that returns a datetime.datetime object.
    Time can be constructed from datetime.datetime: arts.Time(datetime.datetime(...))
    
    Parameters
    ----------
    value : Time
        The Time object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    Time or None
        Modified Time if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Time")
    dialog.setMinimumWidth(400)
    
    layout = QVBoxLayout(dialog)
    
    # Label
    layout.addWidget(QLabel("Select date and time:"))
    
    # DateTime editor widget
    datetime_edit = QDateTimeEdit()
    datetime_edit.setCalendarPopup(True)  # Enable calendar popup for date selection
    datetime_edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss.zzz")  # Show milliseconds
    
    # Convert Time.time (datetime.datetime) to QDateTime
    dt = value.time
    qdt = QDateTime(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, dt.microsecond // 1000)
    datetime_edit.setDateTime(qdt)
    
    layout.addWidget(datetime_edit)
    
    # Info label
    info_label = QLabel("Use arrow keys or mouse to adjust. Calendar popup available for date.")
    info_label.setStyleSheet("color: gray; font-style: italic;")
    layout.addWidget(info_label)
    
    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        # Get the QDateTime and convert back to datetime.datetime
        qdt_result = datetime_edit.dateTime()
        dt_result = datetime.datetime(
            qdt_result.date().year(),
            qdt_result.date().month(),
            qdt_result.date().day(),
            qdt_result.time().hour(),
            qdt_result.time().minute(),
            qdt_result.time().second(),
            qdt_result.time().msec() * 1000  # Convert milliseconds to microseconds
        )
        
        # Create new Time from datetime
        return arts.Time(dt_result)
    
    return None

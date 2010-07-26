import sys,os,glob
from PyQt4 import QtGui,QtCore

class InputWidget(QtGui.QWidget):
  def __init__(self, parent, project):
    QtGui.QWidget.__init__(self,parent)
    self.project = project
    self.projectView = parent

    self.text = QtGui.QTextEdit(self)
    self.text.setContentsMargins(0,0,0,0)
    MyHighlighter(self.text)
    self.setContentsMargins(0,0,0,0)
    hbox = QtGui.QHBoxLayout()
    hbox.setSpacing(2)
    hbox.addWidget(self.text)
    self.setLayout(hbox)
    # Watch project model selection.
    QtCore.QObject.connect(self.project, QtCore.SIGNAL("selectionChanged()"),
        self, QtCore.SLOT("respondToSelection()"))
    self.text.setPlainText(self.project.getInputXMLText())

  @QtCore.pyqtSlot()
  def respondToSelection(self):
    self.text.setPlainText(self.project.getInputXMLText())

class MyHighlighter(QtGui.QSyntaxHighlighter):
  def __init__(self, parent):
    QtGui.QSyntaxHighlighter.__init__(self,parent)
    self.parent = parent
    self.commentFmt = QtGui.QTextCharFormat()
    self.commentFmt.setForeground(QtCore.Qt.gray)
    self.commentFmt.setFontItalic(True)
    #self.commentFmt.setFontWeight(QtGui.QFont.Bold)
    self.IN_COMMENT = 1024

  def highlightBlock(self,text):
    state = self.previousBlockState()
    if state < 0:
      state=0 
    # Finish a comment if we are in one.
    index = 0
    if state & self.IN_COMMENT:
      end = text.indexOf("-->",index)
      if end==-1:
        self.setFormat(index, text.size(), self.commentFmt)
        state = self.IN_COMMENT
      else:
        self.setFormat(index, end-index+3, self.commentFmt)
        index += end+3
        state = 0
    # Find comments.
    index = text.indexOf("<!--",index)
    while index >= 0:
      end = text.indexOf("-->",index+2)
      if end==-1:
        state = self.IN_COMMENT
        self.setFormat(index, text.size()-index, self.commentFmt)
        break
      state = 0
      self.setFormat(index, end-index+3, self.commentFmt)
      index = text.indexOf("<!--",end+2)
    self.setCurrentBlockState(state)

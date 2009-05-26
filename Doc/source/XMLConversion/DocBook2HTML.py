from xml.dom.minidom import parse
from DocBook import DocBookHTMLWriter

document = parse(file('../XML/mmtk-user.xml'))
DocBookHTMLWriter(document, 'MMTK', 'html')

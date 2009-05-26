from xml.dom.minidom import parse
from DocBook import DocBookLatexWriter

document = parse(file('../XML/mmtk-user.xml'))
out = open('mmtk-user.tex', 'w')
DocBookLatexWriter(document, out)
out.close()


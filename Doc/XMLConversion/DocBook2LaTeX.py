from esis_builder import EsisBuilder
from DocBook import DocBookLatexWriter
import sys

esis = EsisBuilder()
esis.feed(open('mmtk-user.esis').read())
document = esis.document

out = open('mmtk-user.tex', 'w')
DocBookLatexWriter(document, out)
out.close()

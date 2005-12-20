from esis_builder import EsisBuilder
from DocBook import DocBookHTMLWriter
import os, sys

esis = EsisBuilder()
esis.feed(open('mmtk-user.esis').read())
document = esis.document

DocBookHTMLWriter(document, 'MMTK', 'html')

# An improved DOM walker. The one in xml.dom assumes that you want to
# iterate over the children of all element nodes, which is not very
# convenient. This one lets you choose between defining start_element
# (and optionally end_element) for automatic iteration and do_element
# for a handler that takes care of all child nodes itself. There are also
# some utility methods for such handlers that locate particular
# subelements and extract data from it.

from xml.dom.minidom import Node
import copy

class DOMWalker:

    def doNode(self, node):
        type = node.nodeType
        if type == Node.ELEMENT_NODE:
            name = node.nodeName

            handler = None
            try:
                handler = getattr(self, "do_"+name)
            except AttributeError:
                pass
            if handler is not None:
                handler(node)
                return

            handler = None
            try:
                handler = getattr(self, "start_"+name)
            except AttributeError:
                pass
            if handler is None:
                self.doElement(node)
            else:
                handler(node)
                for child in node.childNodes:
                    self.doNode(child)
                handler = None
                try:
                    handler = getattr(self, "end_"+name)
                except AttributeError:
                    pass
                if handler is not None:
                    handler(node)
        elif type == Node.TEXT_NODE:
            self.doText(node)
        elif type == Node.COMMENT_NODE:
            self.doComment(node)
        elif type == Node.DOCUMENT_NODE:
            for child in node.childNodes:
                self.doNode(child)
        elif type == Node.DOCUMENT_TYPE_NODE:
            for child in node.childNodes:
                self.doNode(child)
        else:
            self.doOtherNode(node)

    def findChildElements(self, node, element_name):
        elements = []
        for child in node.childNodes:
            if child.nodeType == Node.ELEMENT_NODE and \
               child.nodeName == element_name:
                elements.append(child)
        return elements

    def findTextNodes(self, node):
        nodes = []
        for child in node.childNodes:
            if child.nodeType == Node.TEXT_NODE:
                nodes.append(child)
        return nodes

    def getTextOfElement(self, node, element_name):
        element_name = string.split(element_name, '.')
        while element_name:
            found = 0
            for child in node.childNodes:
                if child.nodeType == Node.ELEMENT_NODE and \
                   child.nodeName == element_name[0]:
                    found = 1
                    break
            if not found:
                return ''
            node = child
            element_name = element_name[1:]
        text_nodes = self.findTextNodes(node)
        text = ''
        for node in text_nodes:
            text = text + node.nodeValue
        return text.encode('L1', 'ignore')

    def doElement(self, node):
        pass

    def doText(self, node):
        pass
    
    def doAttribute(self, node):
        pass

    def doComment(self, node):
        pass

    def doOtherNode(self, node):
        pass


# A simple DOM walker that collects cross-reference labels

class XRefLabels(DOMWalker):

    def __init__(self, document):
        self.labels = {}
        self.doNode(document)

    def doElement(self, node):
        label = node.getAttribute('xreflabel')
        id = node.getAttribute('id')
        if label:
            self.labels[id] = label
        for child in node.childNodes:
            self.doNode(child)

    def __getitem__(self, item):
        try:
            return self.labels[item]
        except KeyError:
            return ""


# The LaTeX writer for DocBook documents. It implements only a subset
# of DocBook, with some extensions that I found necessary for
# documenting Python code. And it produces Python syntax for
# <funcsynopsis>.
#
# Here's the DTD that I use:
#
#  <!DOCTYPE book PUBLIC "-//Norman Walsh//DTD DocBk XML V3.1.4//EN" [
#
#    <!-- Redefinition of paramdef to permit defaultvalue -->
#    <!ELEMENT paramdef (#PCDATA 
#  		  | replaceable 
#  		  | parameter | defaultvalue
#  		  | funcparams)*>
#    <!-- Additional element: defaultvalue -->
#    <!ELEMENT defaultvalue (#PCDATA)>
#    <!-- Additional markup for classes -->
#    <!ELEMENT classdescription (classdef, (para | itemizedlist
#                                           | methoddescription)*)>
#    <!ATTLIST classdescription role CDATA #IMPLIED>
#    <!ENTITY % local.synop.class "|classdescription">
#    <!ELEMENT classdef (#PCDATA | class | classinfo)*>
#    <!ELEMENT class (#PCDATA)>
#    <!ELEMENT classinfo (#PCDATA)>
#    <!ELEMENT methoddescription (methoddef, (void | varargs | paramdef+),
#                                 (para | itemizedlist)*)>
#    <!ELEMENT methoddef (#PCDATA | method)*>
#    <!ELEMENT method (#PCDATA)>
#
#  ]>

import string, sys
from StringIO import StringIO

class Stack:
    def __init__(self):
        self.items = []
    def push(self, item):
        self.items.append(item)
    def pop(self):
        item = self.items[-1]
        del self.items[-1]
        return item
    def top(self):
        return self.items[-1]

class DocBookWriter(DOMWalker):

    def collectText(self, node):
        name = self.file.name
        self.files.push(self.file)
        self.file = StringIO()
        self.file.name = name
        for child in node.childNodes:
            self.doNode(child)
        text = self.file.getvalue()
        self.file = self.files.pop()
        return string.strip(text)


class DocBookLatexWriter(DocBookWriter):

    def __init__(self, document, file):
        self.xreflabels = XRefLabels(document)
        self.file = file
        self.section_level = 0
        self.labelnum = 1
        self.item_name = Stack()
        self.files = Stack()
        self.nclasses = 0

        self.file.write("\\documentclass[12pt]{book}\n\n" +
                        "\\usepackage[latin1]{inputenc}\n" +
                        "\\usepackage{fancyheadings}\n" +
                        "\\parindent=0pt\n" +
                        "\\raggedright\n" +
                        "\\pagestyle{fancyplain}\n" +
                        "\\lhead[\\fancyplain{}{\\sl\\leftmark}]{}\n" +
                        "\\chead{}\n" +
                        "\\rhead[]{\\fancyplain{}{\\sl\\leftmark}}\n" +
                        "\\lfoot[\\rm\\thepage]{}\n" +
                        "\\cfoot{}\n" +
                        "\\rfoot[]{\\rm\\thepage}\n" +
                        "\\addtolength{\\headheight}{3pt}\n" +
                        "\\addtolength{\\footskip}{3pt}\n" +
                        "\\renewcommand{\\chaptermark}[1]{\\markboth{#1}{}}\n"+
                        "\\catcode`\@=11\n" +
                        "\\newcommand\\Tableofcontents{" +
                        "\\chapter*{\\contentsname}\@starttoc{toc}}\n" +
                        "\\usepackage{hyperref}\n" +
                        "\\hypersetup{\n" +
                        "pdfstartview=,\n" +
                        "pdfpagemode=UseOutlines,\n" +
                        #"pdftitle={MMTK User's Guide},\n" +
                        #"pdfauthor={Konrad Hinsen},\n" +
                        #"pdfsubject={User's Guide for MMTK},\n" +
                        #"pdfkeywords={Molecular Modelling,simulation,Python,MMTK}\n" +
                        "}\n\n" +
                        "\\begin{document}\n\n")
        self.doNode(document)
        self.file.write("\\end{document}\n")

    def S(self, text):
        text = string.replace(text, '_', '\\_')
        text = string.replace(text, '^', '\\^')
        text = string.replace(text, '%', '\\%')
        text = string.replace(text, '#', '\\#')
        return text

    def doText(self, node):
        text = node.nodeValue.encode('L1', 'ignore')
        newline_at_end = text[-1] == '\n'
        text = string.split(text, '\n')
        text = filter(lambda s: s, text)
        text = string.join(text, '\n')
        if newline_at_end: text = text + '\n'        
        self.file.write(self.S(text))

    def sectionCommand(self):
        if self.section_level < 3:
            return "\\%ssection" % (self.section_level * "sub")
        elif self.section_level == 3:
            return "\\paragraph"
        elif self.section_level == 4:
            return "\\subparagraph"
        else:
            raise ValueError, "section nesting exceeds LaTeX capabilities"

    def doOtherNode(self, node):
        raise ValueError("Unknown node type %s" % str(node))
        #type = xml.dom.core.NODE_CLASS[node.nodeType].__name__
        #if type != 'ProcessingInstruction':
        #    sys.stderr.write("--> Unknown node %s (type %s)\n"
        #                     % (node.nodeName, type))

    def doElement(self, node):
        sys.stderr.write("--> Unknown element %s\n" % node.nodeName)

    def nop(self, node):
        pass

    start_book = nop
    end_book = nop
    do_title = nop
    start_ulink = nop
    end_ulink = nop

    def do_bookinfo(self, node):
        title = self.findChildElements(node, 'title')[0]
        self.file.write('\\title{%s}\n' % self.collectText(title))
        date = self.findChildElements(node, 'date')[0]
        self.file.write('\\date{%s}\n' % self.collectText(date))
        authors = self.findChildElements(node, 'author')
        text = ''
        for author in authors:
            name = self.getTextOfElement(author, 'firstname') + ' ' + \
                   self.getTextOfElement(author, 'surname')
            text = text + self.S(name) + '\\\\\n'
            affiliations = self.findChildElements(author, 'affiliation')
            for aff in affiliations:
                orgdiv = self.getTextOfElement(aff, 'orgdiv')
                if orgdiv:
                    text = text + self.S(orgdiv) + '\\\\\n'
                orgname = self.getTextOfElement(aff, 'orgname')
                if orgname:
                    text = text + self.S(orgname) + '\\\\\n'
                for a in self.findChildElements(aff, 'address'):
                    text = text + string.strip(self.collectText(a)) + '\n'
                text = text[:-1]
        text = string.split(text, '\n')
        text = filter(lambda s: string.strip(s), text)
        text = string.join(text, '\n')
        self.file.write('\\author{%s}\n\\maketitle\n\n' % text)

    start_street = nop
    def end_street(self, node):
        self.file.write('\\\\ ')
    start_postcode = nop
    def end_postcode(self, node):
        self.file.write(' ')
    start_city = nop
    end_city = nop
    start_country = end_street
    end_country = end_street
    def start_email(self, node):
        self.file.write('E-Mail: ')
    end_email = end_street

    def do_toc(self, node):
        self.file.write('\\tableofcontents\n\n')

    def do_chapter(self, node, title=None):
        if title is None:
            title = self.findChildElements(node, 'title')[0]
            title = self.collectText(title)
        id = node.getAttribute('id')
        if id:
            self.file.write("\\chapter{%s}\n" % title)
            self.file.write("\\hypertarget{%s}{}\n" % id)
            self.file.write('\\label{%s}\n\n' % id)
        else:
            self.file.write("\\chapter{%s}\n\n" % title)
        self.chapter_empty = 1
        for child in node.childNodes:
            self.doNode(child)

    do_preface = do_chapter

##      def do_sect1(self, node):
##          title = self.findChildElements(node, 'title')[0]
##          title = self.collectText(title)
##          id = node.getAttribute('id')
##          if not self.chapter_empty:
##              self.file.write("\\newpage\n")
##          self.chapter_empty = 0
##          self.file.write(self.sectionCommand())
##          if id:
##              self.file.write("*{\\hypertarget{%s}{%s}}\n" % (id, title))
##              self.file.write('\\label{%s}\n' % id)
##          else:
##              self.file.write("*{%s}\n" % title)
##          self.file.write("\n")
##          if title[:7] == 'Module ':
##              self.file.write("\\addcontentsline{toc}{%s}{%s}\n"
##                              % (self.sectionCommand()[1:], title))
##          self.section_level = self.section_level + 1
##          for child in node.childNodes:
##              self.doNode(child)
##          self.section_level = self.section_level - 1
##          self.nclasses = 0

    def do_sect(self, node):
        title = self.findChildElements(node, 'title')[0]
        title = self.collectText(title)
        id = node.getAttribute('id')
        self.file.write(self.sectionCommand())
        if id:
            self.file.write("*{\\hypertarget{%s}{%s}}\n" % (id, title))
            self.file.write('\\label{%s}\n' % id)
        else:
            self.file.write("*{%s}\n" % title)
        self.file.write("\n")
        if self.section_level < 3:
            self.file.write("\\pdfbookmark[%d]{%s}{x%d}\n\n"
                            % ((self.section_level+1),
                               title, self.labelnum))
            self.labelnum = self.labelnum + 1
        self.section_level = self.section_level + 1
        for child in node.childNodes:
            self.doNode(child)
        self.section_level = self.section_level - 1
        self.nclasses = 0
    do_sect2 = do_sect
    do_sect3 = do_sect
    do_sect4 = do_sect
    do_sect5 = do_sect
    do_sect2 = do_sect
    do_simplesect = do_sect

    def do_sect1(self, node):
        if not self.chapter_empty:
            self.file.write("\\newpage\n")
        self.chapter_empty = 0
        self.do_sect(node)

    def start_example(self, node):
        title = self.findChildElements(node, 'title')[0]
        title = self.collectText(title)
        self.file.write(title + '\n\n')

    def start_itemizedlist(self, node):
        self.file.write('\\begin{itemize}\n')
        self.item_name.push('')
    def end_itemizedlist(self, node):
        self.file.write('\\end{itemize}\n')
        self.item_name.pop()


    def start_variablelist(self, node):
        self.file.write('\\begin{description}\n')
        self.item_name.push('')
    def end_variablelist(self, node):
        self.file.write('\\end{description}\n')
        self.item_name.pop()
    start_varlistentry = nop
    def do_term(self, node):
        self.item_name.pop()
        self.item_name.push('[%s]' % self.collectText(node))

    def start_listitem(self, node):
        id = node.getAttribute('id')
        if id:
            self.file.write('\\item\\hypertarget{%s}{%s} '
                            % (id, self.item_name.top()))
        else:
            self.file.write('\\item%s ' % self.item_name.top())

    start_para = nop
    def end_para(self, node):
        self.file.write('\n\n')
    start_simpara = start_para
    end_simpara = end_para

    def start_parameter(self, node):
        self.file.write('{\\sf ')
    def end_parameter(self, node):
        self.file.write('}')

    def start_defaultvalue(self, node):
        self.file.write('={\\tt ')
    def end_defaultvalue(self, node):
        self.file.write('}')

    def start_classdescription(self, node):
        classname = self.S(self.getTextOfElement(node, 'classdef.class'))
        classinfo = self.S(self.getTextOfElement(node, 'classdef.classinfo'))
        title = "Class %s: %s" % (classname, classinfo)
        id = node.getAttribute('id')
        if self.nclasses == 0:
            self.file.write('\\hspace{3mm}\\hrule\\hspace{3mm}\n')
        self.file.write(self.sectionCommand())
        if id:
            self.file.write("*{\\hypertarget{%s}{%s}}\n" % (id, title))
            self.file.write('\\label{%s}\n' % id)
        else:
            self.file.write("*{%s}\n" % title)
        if self.section_level < 3:
            self.file.write("\\pdfbookmark[%d]{%s}{x%d}\n\n"
                            % ((self.section_level+1),
                               title, self.labelnum))
            self.labelnum = self.labelnum + 1
        self.file.write("\n")
        self.nclasses = self.nclasses + 1
        self.methods = 0
    do_classdef = nop
    def end_classdescription(self, node):
        if self.methods != 0:
            self.file.write('\\end{itemize}\n')

    def do_methoddescription(self, node):
        methodname = self.S(self.getTextOfElement(node, 'methoddef.method'))
        text =  '{\\bf\\sf %s}(' % methodname
        parameters = self.findChildElements(node, 'paramdef')
        for p in parameters:
            text = text + self.collectText(p) + ', '
        if text[-2:] == ', ':
            text = text[:-2]
        text = text + ')'
        if self.methods == 0:
            self.file.write('\n\n\medskip\n{\\bf Methods:}\\nopagebreak\n' +
                            '\\nopagebreak\n\\begin{itemize}\n')
        self.methods = self.methods + 1
        self.file.write('\\item ' + text + '\\\\\n')
        for child in node.childNodes:
            self.doNode(child)
    do_methoddef = nop
    do_paramdef = nop
    do_void = nop

    def do_funcsynopsis(self, node):
        id = node.getAttribute('id')
        funcname = self.S(self.getTextOfElement(node, 'funcdef.function'))
        text =  '{\\bf\\sf \\hypertarget{%s}{%s}}(' % (id, funcname)
        parameters = self.findChildElements(node, 'paramdef')
        for p in parameters:
            text = text + self.collectText(p) + ', '
        if text[-2:] == ', ':
            text = text[:-2]
        text = text + ')\\\\\n'
        self.file.write(text)
    
    def start_literal(self, node):
        self.file.write('{\\tt ')
    def end_literal(self, node):
        self.file.write('}')

    def start_emphasis(self, node):
        self.file.write('{\\em ')
    def end_emphasis(self, node):
        self.file.write('}')

    def start_superscript(self, node):
        self.file.write('$^{')
    def end_superscript(self, node):
        self.file.write('}$')

    def start_subscript(self, node):
        self.file.write('$_{')
    def end_subscript(self, node):
        self.file.write('}$')

    def do_programlisting(self, node):
        self.file.write('\\begin{verbatim}\n')
        for text in self.findTextNodes(node):
            text = string.strip(text.nodeValue)
            text = string.split(text, '\n')
            text = filter(lambda s: s, text)
            text = string.join(text, '\n')
            self.file.write(text)
        self.file.write('\n\\end{verbatim}\n')

    def do_xref(self, node):
        linkend = node.getAttribute('linkend')
        self.file.write('\\hyperlink{%s}{%s}'
                        % (linkend, self.S(self.xreflabels[linkend])))
        parts = string.split(linkend, ':')
        if len(parts) > 1:
            if parts[0] == 'Module' or parts[0] == 'Class':
                self.file.write(' (page \\pageref{%s})' % linkend)

    start_bibliography = nop

    def do_bibliodiv(self, node):
        self.nbiblio = 0
        self.do_chapter(node)
        if self.nbiblio > 0:
            self.file.write('\\end{list}\n')

    def do_biblioentry(self, node):
        id = node.getAttribute('id')
        if self.nbiblio == 0:
            self.file.write('\\begin{list}{}{' +
                             '\\renewcommand{\\makelabel}[1]{[#1]}' +
                             '}\n')
        self.nbiblio = self.nbiblio + 1
        self.file.write('\\item[\\hypertarget{%s}{%s}] '
                        % (id, node.getAttribute('xreflabel')))
        authors = self.findChildElements(node, 'author')
        authorgroups = self.findChildElements(node, 'authorgroup')
        for ag in authorgroups:
            authors = authors + self.findChildElements(ag, 'author')
        text = ''
        for author in authors:
            name = self.S(self.getTextOfElement(author, 'firstname')) + ' ' + \
                   self.S(self.getTextOfElement(author, 'surname'))
            text = text + name + ', '
        self.file.write(text[:-2] + '\\\\\n')
        pages = self.findChildElements(node, 'artpagenums')
        if pages:
            title = self.findChildElements(node, 'title')[0]
            text = self.collectText(title) + '\\\\\n'
            text = text + self.S(self.getTextOfElement(node, 'bibliomisc'))
            volume = self.S(self.getTextOfElement(node, 'volumenum'))
            if volume:
                text = text + " {\\bf %s}" % self.S(volume)
            issue = self.getTextOfElement(node, 'issuenum')
            if issue:
                text = text + "(%s)" % self.S(issue)
            text = text + ", %s (%s)\\\\\n" % \
                   (self.collectText(pages[0]),
                    self.S(self.getTextOfElement(node, 'pubdate')))
        else:
            title = self.findChildElements(node, 'title')[0]
            text = self.collectText(title) + '\\\\\n'
            misc = self.findChildElements(node, 'bibliomisc')
            if misc:
                text = text + self.collectText(misc[0]) + "\\\\\n"
            publisher = self.findChildElements(node, 'publisher')
            if publisher:
                text = text + self.S(self.getTextOfElement(publisher[0],
                                                           'publishername')) +\
                              ' ' + \
                              self.S(self.getTextOfElement(node, 'pubdate')) +\
                              '\\\\\n'
            isbn = self.findChildElements(node, 'isbn')
            if isbn:
                text = text + "ISBN: %s\\\\\n" % \
                       self.collectText(isbn[0])
        self.file.write(text)
        self.file.write('\n')

    def do_authorgroup(self, node):
        n = len(node.childNodes)
        for i in range(n):
            self.doNode(node.childNodes[i])
            if i < n-1:
                self.file.write(', ')

    def do_author(self, node):
        text = string.join(map(self.collectText, node.childNodes), ' ')
        self.file.write(text)

    start_firstname = nop
    start_surname = nop
    start_publisher = nop
    start_publishername = nop
    start_pubdate = nop

    def start_isbn(self, node):
        self.file.write('ISBN: ')

    start_bibliomisc = nop
    start_artpagenums = nop
    start_citetitle = nop

    def start_volumenum(self, node):
        self.file.write('{\\bf ')
    def end_volumenum(self, node):
        self.file.write('}')

    def start_issuenum(self, node):
        self.file.write('(')
    def end_issuenum(self, node):
        self.file.write(')')

    def do_glossary(self, node):
        self.do_chapter(node, 'Glossary')
    start_glossentry = nop
    def start_glossterm(self, node):
        id = node.getAttribute('id')
        self.file.write('{\\bf \\hypertarget{%s}{' % id)
    def end_glossterm(self, node):
        self.file.write('}}\\\\\n')
    start_glossdef = nop


# HTML writer for DocBook documents. It is subject to the same
# restrictions as the LaTeX writer, of course.

class HTMLFiles(DOMWalker):

    def __init__(self, document, base):
        self.base_name = base + '_%d.html'
        self.file_number = 0
        self.filenames = {}
        self.files = Stack()
        self.doNode(document)

    def newFileName(self):
        self.file_number = self.file_number + 1
        return self.base_name % self.file_number
    
    def doElement(self, node):
        label = node.getAttribute('xreflabel')
        id = node.getAttribute('id')
        if label:
            self.filenames[id] = self.files.top()
        for child in node.childNodes:
            self.doNode(child)

    def do_chapter(self, node):
        self.files.push(self.newFileName())
        node.setAttribute('filename', self.files.top())
        self.doElement(node)

    do_bibliodiv = do_chapter
    do_glossary = do_chapter

    def do_sect(self, node):
        title = self.getTextOfElement(node, 'title')
        self.doElement(node)

    do_preface = do_sect
    do_sect1 = do_sect
    do_sect2 = do_sect
    do_sect3 = do_sect
    do_sect4 = do_sect
    do_sect5 = do_sect
    do_simplesect = do_sect

    def __getitem__(self, item):
        try:
            return self.filenames[item]
        except KeyError:
            return ""


class DocBookHTMLWriter(DocBookWriter):

    def __init__(self, document, filename, directory=''):
        self.xreflabels = XRefLabels(document)
        if directory:
            self.directory = directory+'/'
        else:
            self.directory = directory
        self.filenames = HTMLFiles(document, filename)
        self.file = open(self.directory+filename+'.html', 'w')
        self.files = Stack()
        self.section_levels = Stack()
        self.section_levels.push(3)
        self.submodules = {}
        self.doNode(document)
        self.file.close()

    def doText(self, node):
        text = node.nodeValue
        initial_newline = text and text[0] == '\n'
        terminal_newline = text and text[-1] == '\n'
        text = string.split(text, '\n')
        text = filter(lambda s: s, text)
        text = string.join(text, '\n')
        if initial_newline: text = '\n' + text
        if terminal_newline: text = text + '\n'        
        text = text.encode('L1', 'ignore')
        self.file.write(text)

    def doOtherNode(self, node):
        type = xml.dom.core.NODE_CLASS[node.nodeType].__name__
        if type != 'ProcessingInstruction':
            sys.stderr.write("--> Unknown node %s (type %s)\n"
                             % (node.nodeName, type))

    def doElement(self, node):
        sys.stderr.write("--> Unknown element %s\n" % node.nodeName)

    def nop(self, node):
        pass

    start_book = nop
    end_book = nop
    do_title = nop
    do_toc = nop

    def do_bookinfo(self, node):
        title = self.findChildElements(node, 'title')[0]
        title = self.collectText(title)
        self.file.write('<center>\n')
        self.file.write('<title>%s</title>\n' % title)
        self.file.write('<h1>%s</h1>\n' % title)
        authors = self.findChildElements(node, 'author')
        text = ''
        for author in authors:
            name = self.getTextOfElement(author, 'firstname') + ' ' + \
                   self.getTextOfElement(author, 'surname')
            text = text + name + '<br>\n'
            affiliations = self.findChildElements(author, 'affiliation')
            for aff in affiliations:
                orgdiv = self.getTextOfElement(aff, 'orgdiv')
                if orgdiv:
                    text = text + orgdiv + '<br>\n'
                orgname = self.getTextOfElement(aff, 'orgname')
                if orgname:
                    text = text + orgname + '<br>\n'
                for a in self.findChildElements(aff, 'address'):
                    text = text + string.strip(self.collectText(a)) + '\n'
                text = text[:-1]
        text = string.split(text, '\n')
        text = filter(lambda s: string.strip(s), text)
        text = string.join(text, '\n')
        self.file.write(text)
        date = self.findChildElements(node, 'date')[0]
        self.file.write('Last revision: %s<br>' % self.collectText(date))
        self.file.write('</center>\n')

    def newline(self, node):
        self.file.write('<br>\n')
        
    start_street = nop
    end_street = newline
    start_postcode = nop
    def end_postcode(self, node):
        self.file.write(' ')
    start_city = nop
    end_city = nop
    start_country = newline
    end_country = newline

    def do_email(self, node):
        email = self.collectText(node)
        self.file.write('E-Mail: <a href="mailto:%s">%s</a><br>\n'
                        % (email, email))

    def do_sect(self, node, chapter=0, title=None):
        if not title:
            title = self.findChildElements(node, 'title')[0]
            title = self.collectText(title)
        id = node.getAttribute('id')
        if chapter or title[:7] == 'Module ':
            filename = node.getAttribute('filename')
            if not chapter and not self.submodules.get(self.file, 0):
                self.file.write('<h1>Submodules:</h1>\n')
                self.submodules[self.file] = 1
            self.file.write('<ul><li><a href="%s">%s</a></ul>\n'
                            % (filename, title))
            self.files.push(self.file)
            self.file = open(self.directory+filename, 'w')
            self.section_levels.push(2)
            if id:
                self.file.write('<a name="%s"><h1>%s</h1></a>\n\n'%(id,title))
            else:
                self.file.write("<h1>%s</h1>\n\n" % title)
            children = list(copy.copy(node.childNodes))
            local = []
            for i in range(len(children)):
                child = children[i]
                if child.nodeType == Node.ELEMENT_NODE:
                    if child.nodeName[:4] == 'sect':
                        text = self.collectText(child)
                        if text[:9] == 'Functions':
                            local.append(child)
                            children[i] = None
                    elif child.nodeName == 'classdescription':
                        local.append(child)
                        children[i] = None
                if children[i] is not None:
                    self.doNode(child)
            for child in local:
                self.doNode(child)
            self.file.close()
            self.section_levels.pop()
            self.file = self.files.pop()
        else:
            if self.section_levels.top() == 2:
                self.file.write('<hr width=70%>\n')
            if id:
                self.file.write('<a name="%s"><h%d>%s</h%d><a>\n\n' %
                                (id, self.section_levels.top(),
                                 title, self.section_levels.top()))
            else:
                self.file.write("<h%d>%s</h%d>\n\n" %
                                (self.section_levels.top(),
                                 title, self.section_levels.top()))
            self.section_levels.push(self.section_levels.pop() + 1)
            for child in node.childNodes:
                self.doNode(child)
            self.section_levels.push(self.section_levels.pop() - 1)

    do_preface = do_sect
    do_sect1 = do_sect
    do_sect2 = do_sect
    do_sect3 = do_sect
    do_sect4 = do_sect
    do_sect5 = do_sect
    do_simplesect = do_sect

    def do_chapter(self, node):
        self.do_sect(node, 1)

    def start_example(self, node):
        title = self.findChildElements(node, 'title')[0]
        title = self.collectText(title)
        self.file.write(title + '\n\n')

    def start_itemizedlist(self, node):
        self.file.write('<ul>\n')
        self.listtype = 'unordered'
    def end_itemizedlist(self, node):
        self.file.write('</ul>\n')


    def start_variablelist(self, node):
        self.file.write('<dl>\n')
        self.listtype = 'definition'
    def end_variablelist(self, node):
        self.file.write('</dl>\n')
    start_varlistentry = nop
    def do_term(self, node):
        self.file.write('<dt>%s</dt>\n' % self.collectText(node))

    def do_listitem(self, node):
        if self.listtype == 'definition':
            self.file.write('<dd>%s</dd>\n' % self.collectText(node))
        else:
            id = node.getAttribute('id')
            if id:
                self.file.write('<li><a name="%s" %s</a>' %
                                (id, self.collectText(node)))
            else:
                self.file.write('<li> %s' % self.collectText(node))

    def start_para(self, node):
        self.file.write('<p>')
    def end_para(self, node):
        self.file.write('</p>\n\n')
    start_simpara = start_para
    end_simpara = end_para

    def start_parameter(self, node):
        self.file.write('<i>')
    def end_parameter(self, node):
        self.file.write('</i>')

    def start_filename(self, node):
        self.file.write('<i>')
    def end_filename(self, node):
        self.file.write('</i>')

    def start_defaultvalue(self, node):
        self.file.write('=<tt>')
    def end_defaultvalue(self, node):
        self.file.write('</tt>')

    def start_classdescription(self, node):
        classname = self.getTextOfElement(node, 'classdef.class')
        classinfo = self.getTextOfElement(node, 'classdef.classinfo')
        id = node.getAttribute('id')
        title = "Class %s: %s" % (classname, classinfo)
        self.file.write('<hr width=70%>\n')
        if id:
            self.file.write('<a name="%s"><h%d>%s</h%d></a>\n\n'
                            % (id, self.section_levels.top(),
                               title, self.section_levels.top()))
        else:
            self.file.write('<h%d>%s</h%d>\n\n'
                            % (self.section_levels.top(), title,
                               self.section_levels.top()))
        self.methods = 0
    do_classdef = nop
    def end_classdescription(self, node):
        if self.methods != 0:
            self.file.write('</ul>\n')

    def do_methoddescription(self, node):
        methodname = self.getTextOfElement(node, 'methoddef.method')
        text =  '<b><i>%s</i></b>(' % methodname
        parameters = self.findChildElements(node, 'paramdef')
        for p in parameters:
            text = text + self.collectText(p) + ', '
        if text[-2:] == ', ':
            text = text[:-2]
        text = text + ')'
        if self.methods == 0:
            self.file.write('<b>Methods:</b><br>\n<ul>\n')
        self.methods = self.methods + 1
        self.file.write('<li> ' + text + '\n')
        for child in node.childNodes:
            self.doNode(child)
    do_methoddef = nop
    do_paramdef = nop
    do_void = nop

    def do_funcsynopsis(self, node):
        funcname = self.getTextOfElement(node, 'funcdef.function')
        id = node.getAttribute('id')
        text = ''
        if id:
            text =  '<a name="%s">' % id
        text =  text + '<b><i>%s</i></b>(' % funcname
        parameters = self.findChildElements(node, 'paramdef')
        for p in parameters:
            text = text + self.collectText(p) + ', '
        if text[-2:] == ', ':
            text = text[:-2]
        text = text + ')'
        if id:
            text = text + '</a>'
        text = text + '<br>\n'
        self.file.write(text)

    def start_literal(self, node):
        self.file.write('<tt>')
    def end_literal(self, node):
        self.file.write('</tt>')

    def start_emphasis(self, node):
        self.file.write('<i>')
    def end_emphasis(self, node):
        self.file.write('</i>')

    def start_superscript(self, node):
        self.file.write('<sup>')
    def end_superscript(self, node):
        self.file.write('</sup>')

    def start_subscript(self, node):
        self.file.write('<sub>')
    def end_subscript(self, node):
        self.file.write('</sub>')

    def do_programlisting(self, node):
        self.file.write('<pre>\n')
        for text in self.findTextNodes(node):
            text = string.strip(text.nodeValue)
            text = string.split(text, '\n')
            text = filter(lambda s: s, text)
            text = string.join(text, '\n')
            self.file.write(text)
        self.file.write('\n</pre>\n')

    def start_ulink(self, node):
        self.file.write('<a href="%s">' % node.getAttribute('url'))
    def end_ulink(self, node):
        self.file.write('</a>')

    def do_xref(self, node):
        linkend = node.getAttribute('linkend')
        filename = self.filenames[linkend]
        if filename == self.file.name:
            filename = ''
        self.file.write('<a href="%s#%s">' % (filename, linkend))
        self.file.write(self.xreflabels[linkend])
        self.file.write('</a>')

    start_bibliography = nop

    def do_bibliodiv(self, node):
        self.nbiblio = 0
        self.do_chapter(node)
        if self.nbiblio > 0:
            self.file.write('</ul>\n')

    def do_biblioentry(self, node):
        if self.nbiblio == 0:
            self.file.write('<ul>\n')
        self.nbiblio = self.nbiblio + 1
        self.file.write('<li><b>[%s]</b><br>\n' % node.getAttribute('xreflabel'))
        authors = self.findChildElements(node, 'author')
        authorgroups = self.findChildElements(node, 'authorgroup')
        for ag in authorgroups:
            authors = authors + self.findChildElements(ag, 'author')
        text = ''
        for author in authors:
            name = self.getTextOfElement(author, 'firstname') + ' ' + \
                   self.getTextOfElement(author, 'surname')
            text = text + name + ', '
        self.file.write(text[:-2] + '<br>\n')
        pages = self.findChildElements(node, 'artpagenums')
        if pages:
            title = self.findChildElements(node, 'title')[0]
            text = self.collectText(title) + '<br>\n'
            text = text + self.getTextOfElement(node, 'bibliomisc')
            volume = self.getTextOfElement(node, 'volumenum')
            if volume:
                text = text + " <b>%s</b>" % volume
            issue = self.getTextOfElement(node, 'issuenum')
            if issue:
                text = text + "(%s)" % issue
            text = text + ", %s (%s)<br>\n" % (self.collectText(pages[0]),
                              self.getTextOfElement(node, 'pubdate'))
        else:
            title = self.findChildElements(node, 'title')[0]
            text = self.collectText(title) + '<br>\n'
            misc = self.findChildElements(node, 'bibliomisc')
            if misc:
                text = text + self.collectText(misc[0]) + "<br>\n"
            publisher = self.findChildElements(node, 'publisher')
            if publisher:
                text = text + self.getTextOfElement(publisher[0],
                                                    'publishername') + ' ' + \
                              self.getTextOfElement(node, 'pubdate') + '<br>\n'
            isbn = self.findChildElements(node, 'isbn')
            if isbn:
                text = text + "ISBN: %s<br>\n" % self.collectText(isbn[0])
        self.file.write(text)

##      def do_biblioentry(self, node):
##          if self.nbiblio == 0:
##              self.file.write('<ol>\n')
##          self.nbiblio = self.nbiblio + 1
##          self.file.write('<li>')
##          authors = self.findChildElements(node, 'author')
##          authorgroups = self.findChildElements(node, 'authorgroup')
##          for ag in authorgroups:
##              authors = authors + self.findChildElements(ag, 'author')
##          text = ''
##          for author in authors:
##              name = self.getTextOfElement(author, 'firstname') + ' ' + \
##                     self.getTextOfElement(author, 'surname')
##              text = text + name + ', '
##          self.file.write(text[:-2] + '<br>\n')
##          pages = self.findChildElements(node, 'artpagenums')
##          if pages:
##              title = self.findChildElements(node, 'title')[0]
##              text = self.collectText(title) + '<br>\n'
##              text = text + self.getTextOfElement(node, 'bibliomisc')
##              volume = self.getTextOfElement(node, 'volumenum')
##              if volume:
##                  text = text + " <b>%s</b>" % volume
##              issue = self.getTextOfElement(node, 'issuenum')
##              if issue:
##                  text = text + "(%s)" % issue
##              text = text + ", %s (%s)<br>\n" % (self.collectText(pages[0]),
##                                self.getTextOfElement(node, 'pubdate'))
##          else:
##              title = self.findChildElements(node, 'title')[0]
##              text = self.collectText(title) + '<br>\n'
##              misc = self.findChildElements(node, 'bibliomisc')
##              if misc:
##                  text = text + self.collectText(misc[0]) + "<br>\n"
##              publisher = self.findChildElements(node, 'publisher')
##              if publisher:
##                  text = text + self.getTextOfElement(publisher[0],
##                                                      'publishername') + ' ' + \
##                                self.getTextOfElement(node, 'pubdate') + '<br>\n'
##              isbn = self.findChildElements(node, 'isbn')
##              if isbn:
##                  text = text + "ISBN: %s<br>\n" % self.collectText(isbn[0])
##          self.file.write(text)

    def do_authorgroup(self, node):
        n = len(node.childNodes)
        for i in range(n):
            self.doNode(node.childNodes[i])
            if i < n-1:
                self.file.write(', ')

    def do_author(self, node):
        text = string.join(map(self.collectText, node.childNodes), ' ')
        self.file.write(text)

    start_firstname = nop
    start_surname = nop
    start_publisher = nop
    start_publishername = nop
    start_pubdate = nop

    def start_isbn(self, node):
        self.file.write('ISBN: ')

    start_bibliomisc = nop
    start_artpagenums = nop
    start_citetitle = nop

    def start_volumenum(self, node):
        self.file.write('<b>')
    def end_volumenum(self, node):
        self.file.write('</b>')

    def start_issuenum(self, node):
        self.file.write('(')
    def end_issuenum(self, node):
        self.file.write(')')

    def do_glossary(self, node):
        self.do_sect(node, 1, 'Glossary')
    start_glossentry = nop
    def start_glossterm(self, node):
        self.file.write('<a name="%s"><b>' % node.getAttribute('id'))
    def end_glossterm(self, node):
        self.file.write('</b></a><br>\n')
    start_glossdef = nop

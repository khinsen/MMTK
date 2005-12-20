'''DOM builder from ESIS output.'''

# Hacked from ``pygrove.py'' from Paul Prescod.

import string, re, regsub, sys
from xml.dom.core import *
from xml.dom.builder import Builder


_sdata_dict = {
    'acirc': 'â',
    'agrave': 'à',
    'ccedil': 'ç',
    'dollar': '$',
    'eacute': 'é',
    'ecirc': 'ê',
    'egrave': 'è',
    'icirc': 'î',
    'ocirc': 'ô',
    'percnt': '%',
    'ugrave': 'ù',
}

def handle_sdata(sdata):
    return _sdata_dict.get(sdata, 'unknown')


class EsisBuilder(Builder):
    
    def __init__(self):
        Builder.__init__(self)
        self.attr_store = {}
        self.id_store = {}
        #self.sdata_handler = handle_sdata

    def feed(self, data):
        for line in string.split(data, '\n'):
            if not line: 
                break
            event = line[0]
            text = line[1:]

            if event == '(':
                element = self.document.createElement(text, self.attr_store)
                self.attr_store = {}
                self.push(element)

            elif event == ')':
                self.pop()

            elif event == 'A':
                l = re.split(' ', text, 2)
                if len(l) > 2:
                    name = l[0]
                    value = ESISDecode(l[2])
                    self.attr_store[name] = value

            elif event == '-':
                text = self.document.createText(ESISDecode(text))
                self.push(text)
            elif event == '?':
                text = string.split( text[:-1] )
                pi = self.document.createProcessingInstruction(text[0],
                                           string.join(text[1:]) )
                self.push( pi )
                
            elif event == 'C':
                return

            elif event == 'e':
                # Indicates that this is an empty element;
                # only produced by nsgmls for -oempty.    We
                # can safely ignore it.
                pass

            elif event in 'spf':
                # Some sort of command that applies to a
                # following command; save it 
                self.id_store[ event ] = text

            elif event == 'N':
                pubId = sysId = ""
                if self.id_store.has_key('p'): 
                    pubId = self.id_store['p']
                if self.id_store.has_key('s'): 
                    sysId = self.id_store['s']
                notation = self.document.createNotation(text,
                                    pubId,
                                    sysId)
                self.id_store = {}
            
            else:
                sys.stderr.write('Unknown event: ' + `line` + '\n')


backslash = r"\\"
regor = "|"

find = "(" + r"\\\\" + regor \
        + r"\\n" + regor \
        + r"\\\|\[[^" +backslash+ "]*" + r"\]\\\|" + regor \
        + r"\\[0-9]+" + regor \
        + r"\\#[0-9]+" + regor \
        + r"\\%[0-9]+" + regor \
        + ")"

def fix(text):
    if (not text) or (text[0] != "\\"):
        return text
    text = text[1:]
    if (text == "\\"):
        return "\\"
    elif (text == "n"):
        return "\n"
    elif(text[0] == "|"):
        return handle_sdata(string.strip(text[2:-3]))
        #return '&' + string.strip(text[2:-3]) + ';'
    elif(text[0] == "#"):
        return chr(string.atoi(text[1:],10))
    elif(text[0] == "%"):
        return chr(string.atoi(text[1:],10))
    else:
        return chr(string.atoi(text[0:3], 8)) + text[3:]

def ESISEncode(text):
    return re.sub("\n", "\\n", text);

def ESISDecode(text):
    #return regsub.gsub("\\\\n", "\n", text)

    #print `text`
    prog = re.compile(find)
    parts = prog.split(text, find)
    #print parts
    parts = map(fix, parts)
    res = string.join(parts, "");
    #print res
    return res


if __name__ == '__main__':
    import sys
    from xml.dom.writer import XmlLineariser

    p = EsisBuilder()
    p.feed(open(sys.argv[1]).read())
    w = XmlLineariser()
    w.add_newline_after = [ 'p', 'title', 'abstract' ]
    print w.linearise(p.document.documentElement)

# vim:ts=2:ai

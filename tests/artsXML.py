import xml.sax
from numpy import *
#from general import *
import os
import cStringIO
import types
from warnings import warn

arts_dim_labels=["ncols","nrows","npages","nbooks","nshelves","nvitrines","nlibraries"]

arts_tensor_names=["Vector","Matrix","Tensor3","Tensor4","Tensor5","Tensor6","Tensor7"]

arts_text_names=["String","comment","SpeciesTag"]

verbosity=0#controls the verbosity of all classes and functions in this module

class Tensor:
    def __init__(self,rank,attributes):
	if verbosity:
	    print "creating Tensor of rank " + str(rank)
	self.dimlist = []
	if (rank==1):
	    self.dimlist.append(int(attributes.getValue("nelem")))
	else:
	    for i in range(rank):
		self.dimlist.append(int(attributes.getValue(arts_dim_labels[i])))
	#put the highest dimensions first	
	self.dimlist.reverse()
	if verbosity:
	    print "dimensions:"
	    print self.dimlist
    def fixshape(self):
	if verbosity:
	    print "raw data shape = "
	    print self.array.shape
        self.array.shape=self.dimlist

class newHandler(xml.sax.ContentHandler):
    """XML handler for arts files.  This is much more general then the old one.
    Unlike the previous handler, this handler works by continutally
    updating a list of strings that constitute a reference for the eventual
    location of an object within the overall data structure. With the exception
    of Arrays, which are represented by lists, every XML tag becomes a
    dictionary key.  In cases where the tag line has a name attribute, the
    dictionary key becomes the name string. This behaviour is enabled by default
    but can be controlled with the une_names argument on initialisation.

    Normally it should not be necessary to use this function directly, instead
    use the newLoad function below.
    """
    def __init__(self,use_names=True):
        """If use_names is True, the handler will use any present \"name\"
        attributes as keys for each object. Otherwise tags are used"""
        self.CurrentHierarchy=[]
        self.text=""
        self.datastruct={}
        self.reference_str_list=["self.datastruct"]
	self.assignment_str_list=[]#this is used to ensure nothing is overwritten
        self.use_names=use_names
    def characters(self,text):
        if self.StringIOObj_isopen:
            self.StringIOobj.write(text)
    def startDocument(self):
	if verbosity:
	    print "Starting to parse XML"
    def startElement(self,tag,attributes):
	def advance_name_if_nec(name,reference_str_list):
	    if (eval(reduce(lambda x, y: x+y, self.reference_str_list)+\
                 ".has_key(\""+name+"\")")):
		#then change the name
		if verbosity:
		    print 'already have '+name
		if (len(name.split())==1):
		    name=name+' 0'
		else:
		    label=name.split()[0]
		    number=int(name.split()[1])
		    name=label+' '+str(number+1)
		if verbosity:
		    print "name changed to "+name
		name=advance_name_if_nec(name,reference_str_list)
	    return name
	self.StringIOobj=cStringIO.StringIO()
        self.StringIOObj_isopen=1
        if self.use_names:
            name=attributes.get("name",tag)
        else:
            name=tag
	if verbosity:
	    print "begin Name: "+name+"  Tag: "+tag
        #Use tag to modify reference string for subsequent data
        if tag=="Array":
            name=advance_name_if_nec(name,self.reference_str_list)
	    self.CurrentHierarchy.append(name)
        
	    list_of_empty_dicts=[];#[{}]*n gives n references to the same dict!!!!
	    for i in range(int(attributes["nelem"])):
		list_of_empty_dicts.append({})
	    eval(reduce(lambda x, y: x+y, self.reference_str_list)+\
                 ".update({\""+name+"\":list_of_empty_dicts})")
            self.reference_str_list.append("[\""+name+"\"]")
            self.reference_str_list.append("[0]")
            self.CurrentHierarchy.append(0)
        else:
	    #check that new object doesn't already exist
	    
	    name=advance_name_if_nec(name,self.reference_str_list)
	    
            eval(reduce(lambda x, y: x+y, self.reference_str_list)+\
                 ".update({\""+name+"\":{}})")
            self.reference_str_list.append("[\""+name+"\"]")
	    self.CurrentHierarchy.append(name)
        
	    #If appropriate, create a tensor object
        if tag in arts_tensor_names:
            rank = arts_tensor_names.index(tag)+1
            self.CurrentObject=Tensor(rank,attributes)
    def endElement(self,tag):
        if self.StringIOObj_isopen:
            self.text=self.StringIOobj.getvalue()
            self.StringIOobj.close()
            self.StringIOObj_isopen=0
	if verbosity:
	    print "end  Tag: "+tag
        if len(self.text.strip())>0:
            if tag in arts_tensor_names:
                splittext=self.text.split()
                datalist=[]
                if len(splittext)>0:
                    for i in range(len(splittext)):
                        datalist.append(float(splittext[i]))
                self.CurrentObject.array=array(datalist)
                self.CurrentObject.fixshape()
                self.CurrentObject=self.CurrentObject.array
            elif tag in arts_text_names:
                self.CurrentObject=self.text.strip()
            elif tag=="Index":
                self.CurrentObject=int(self.text.strip())
            else:
                self.CurrentObject=float(self.text.strip())
	    assignment_str=reduce(lambda x, y: x+y,
                        self.reference_str_list[:len(self.reference_str_list)-1])+\
                 ".update({\""+self.CurrentHierarchy[len(self.CurrentHierarchy)-1]+\
                 "\":self.CurrentObject})"
	    if verbosity:
		print assignment_str
            eval(assignment_str)
        self.text=""
        
        self.CurrentHierarchy.pop()
        self.reference_str_list.pop()
        if tag=="Array":
            self.CurrentHierarchy.pop()
            self.reference_str_list.pop()
	#if necessary advance array index
	if len(self.CurrentHierarchy)>0:
	    ref=self.CurrentHierarchy[len(self.CurrentHierarchy)-1]
	    if (type(ref)==types.IntType):
		ref+=1
		self.CurrentHierarchy[len(self.CurrentHierarchy)-1]=ref
		self.reference_str_list[len(self.CurrentHierarchy)]="["+str(ref)+"]"


def load(filename,use_names=True):
    """This general purpose function returns a dictionary structure reflecting
    the structure of the XML file.  If there is only one object in the structure,
    then that single object is returned. As far as far as I know, this works with
    every data type exported by ARTS

    Usage:
          data_struct=load(filename)

    """
    if not os.path.exists(filename):
        raise IOError,'file: '+filename+' not found'
    import xml.sax
    handler=newHandler(use_names)
    xml.sax.parse(filename,handler)
    data=handler.datastruct["arts"]
    if len(data.keys())==1:
        return data.values()[0]
    else:
        return data



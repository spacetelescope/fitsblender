""" blendheaders - Merge headers from multiple inputs to create a new header and table

"""
import os
import numpy as np
import pyfits

from stsci.tools import fileutil, textutil, parseinput

from . import blender

__taskname__ = 'blendheaders' # unless someone comes up with anything better
__version__ = '0.1.1'
__vdate__ = '22-Feb-2012'

fits_required_bool_kws = ['SIMPLE','EXTEND']
WCS_KEYWORDS=['CD1_1','CD1_2', 'CD2_1', 'CD2_2', 'CRPIX1',
'CRPIX2','CRVAL1', 'CRVAL2', 'CTYPE1', 'CTYPE2',
'VAFACTOR','ORIENTAT','BUNIT','WCSNAME']

#### Custom blending functions
def multi(vals):
    """
    This will either return the common value from a list of identical values
    or 'MULTIPLE' 
    """
    uniq_vals = list(set(vals))
    num_vals = len(uniq_vals)
    if num_vals == 0:
        return None
    if num_vals == 1:
        return uniq_vals[0]
    if num_vals > 1:
        return "MULTIPLE"

# translation dictionary for function entries from rules files
blender_funcs = {'first':blender.first,
                'last':blender.last,
                'multi':multi,
                'mean':np.mean,
                'sum':np.sum,
                'max':np.max,
                'min':np.min,
                'stddev':np.std}

delete_command = '<delete>'

#### TEAL Interfaces to run this task    

def getHelpAsString(docstring=False):
    """ 
    return useful help from a file in the script directory called __taskname__.help
    """
    install_dir = os.path.dirname(__file__)
    htmlfile = os.path.join(install_dir,'htmlhelp',__taskname__+'.html')
    helpfile = os.path.join(install_dir,__taskname__+'.help')
    if docstring or (not docstring and not os.path.exists(htmlfile)):
        helpString = __taskname__+' Version '+__version__+' updated on '+__vdate__+'\n\n'
        if os.path.exists(helpfile):
            helpString += teal.getHelpFileAsString(__taskname__,__file__)
        else:
            helpString += blendheaders.__doc__
    else:
        helpString = 'file://'+htmlfile

    return helpString

def run(configobj):
    # Interpret user-input from TEAL GUI and call function
    blendheaders(configobj['drzfile'],
                inputs=configobj['inputs'],
                output=configobj['output'],
                sciext=configobj['sciext'],
                errext=configobj['errext'],
                dqext=configobj['dqext'],
                verbose=configobj['verbose'])

#### Primary functional interface for the code
def blendheaders(drzfile, inputs=None, output=None, 
                sciext='SCI', errext='ERR', dqext='DQ', 
                verbose=False):
    """ Blend headers that went into creating the original drzfile into a 
    new header with table that contains keyword values from all input images.
    
    The drzfile will be used to determine the names of the input files, should
    no filenames be provided in the 'inputs' parameter.
    
    The drzfile will be updated 'in-place' with the new headers and table if
    no output filename has been provided.
    
    Parameters
    ----------
    drzfile : str
        Name of drizzled image(s) with headers that need updating. This can
        be specified as a single filename, or using wildcards, or '@'-file or
        python list of filenames.
        When no value for 'inputs' has been provided, this file(or set of files) 
        will be used to determine the names of the input (flt.fits) files 
        whose headers need to be blended to create the new drzfile header.
    
    inputs : list, optional
        If provided, the filenames with extensions for each chip 
        provided in this list will be used to get
        the headers which will be blended into the final output headers.
        For example, ['j9cd01kqq_flt.fits[sci,1]','j9cd01kqq_flt.fits[sci,2]']
        would create a blended header based on these chips headers.

    output : str, optional
        If specified, a new file will be written out that contains the updated
        (blended) headers.
    
    sciext: str, optional [Default: 'SCI']
        EXTNAME of extensions with science data from the input FITS files. The 
        header of this extension will be used as the basis for the SCI header
        of the drizzled product FITS file.

    errext: str, optional [Default: 'ERR']
        EXTNAME of extensions with the error array from the input FITS files. The 
        header of this extension will be used as the basis for the WHT header
        of the drizzled product FITS file. If blank or "INDEF", it will use
        the 'SCI' header as the basis for the output header for the WHT array.

    dqext: str, optional [Default: 'DQ']
        EXTNAME of extensions with the data quality array from the input FITS 
        files. The header of this extension will be used as the basis for the 
        CTX header of the drizzled product FITS file (when a CTX extension gets
        created). If blank or "INDEF", it will use the 'SCI' header as the 
        basis for the header of any generated CTX array.
        
    verbose : bool, optional [Default: False]
        Print out additional messages during processing when specified.
        
    """
    # interpret input 
    drzfiles = parseinput.parseinput(drzfile)[0]

    # operate on each drzfile specified
    for drzfile in drzfiles:
        if inputs in [None, '',' ','INDEF','None']:
            inputs = extract_filenames_from_drz(drzfile)

        if verbose:
            print 'Creating blended headers from: '
            for i in inputs: print '    ',i

        newhdrs, newtab = get_blended_headers(inputs,verbose=verbose)

                    
        # Remove distortion related keywords not included in rules
        for hdr in newhdrs:
            remove_distortion_keywords(hdr)
        
        # open drizzle product to update headers with new headers
        open_mode='update'
        if output not in [None,'',' ','INDEF','None']:
            open_mode = 'readonly'
        drzimg = pyfits.open(drzfile,mode=open_mode)
        
        # Determine whether we are working with a simple DRZ FITS file or
        # a full multi-extension DRZ FITS file.
        if len(drzimg) < 3:
            # We are working with a simple FITS image, so concatenate the 
            # blended PRIMARY and SCI  headers
            drzimg[0].header = cat_headers(newhdrs[0],newhdrs[1])
            drzimg.append(newtab)
        else:
            # We are working with a full MEF file, so update all extension headers
            for i,extn in enumerate(drzimg):
                if isinstance(extn, pyfits.BinTableHDU):
                    break
                # Update new headers with correct array sizes
                if isinstance(extn, pyfits.ImageHDU):
                    extn_naxis = extn.header['NAXIS']
                    newhdrs[i].update('NAXIS',extn_naxis)
                    newhdrs[i].update('BITPIX',extn.header['BITPIX'])
                    for card in newhdrs[i]['naxis*']:
                        if len(card.key) > 5: # naxisj keywords
                            if extn_naxis > 0: 
                                newhdrs[i].update(card.key, extn.header[card.key])
                            else:
                                del newhdrs[i][card.key]
                    newhdrs[i].update('EXTNAME', extn.header['EXTNAME'],after='ORIGIN')
                    newhdrs[i].update('EXTVER', extn.header['EXTVER'],after='EXTNAME')
                    for kw in WCS_KEYWORDS:
                        if kw in extn.header:
                            newhdrs[i].update(kw,extn.header[kw])
                if isinstance(extn, pyfits.PrimaryHDU):
                    for card in extn.header['exp*']:
                        newhdrs[i].update(card.key,card.value)
                    newhdrs[i].update('NEXTEND',len(drzimg)-1)
                    newhdrs[i].update('ROOTNAME',extn.header['rootname'])
                    newhdrs[i].update('BITPIX',extn.header['bitpix'])
                    # Determine which keywords are included in the table but not 
                    # the new dict(header). These will be removed from the output
                    # header altogether
                    tabcols = newtab.data.dtype.names
                    hdrkws = newhdrs[i].keys()
                    del_kws = list(set(tabcols) - set(hdrkws))
                    del_kws.append('HISTORY')
                    for kw in extn.header:
                        if kw not in newhdrs[i] and kw not in del_kws:
                            newhdrs[i].update(kw,extn.header[kw])

                extn.header = newhdrs[i]

        # Now append table with remaining header keyword values
        drzimg.append(newtab)

        # Write out the updated product
        if open_mode == 'update':
            drzimg.close()
            print 'Updated ',drzfile,' with blended headers.'
        else:
            if os.path.exists(output): os.remove(output)
            drzimg.writeto(output)
            drzimg.close()
            print 'Created new file ',output,' with blended headers.'

        # Clean up for the next run
        del drzimg, newhdrs, newtab


def get_blended_headers(inputs, verbose=False,extlist=['SCI','ERR','DQ']):
    """ 
    Return a set of blended headers based on the input files/headers provided

    Parameters
    ----------
    inputs : list
        Either a single list of filenames from which to extract the headers to
        be blended, or a list of lists of pyfits.Header objects to be blended.
        
    Returns
    -------
    headers : list
        A list of blended headers, one for each image extension recognized in
        the input image list or one for each set of headers provided as input

    new_table : object
        Single pyfits.TableHDU object that contains the combined results from
        all input headers(extension). Each row will correspond to an image, 
        and each column corresponds to a single keyword listed in the rules.

    """    
    hdrlist = None
    if not isinstance(inputs, list):
        inputs = [inputs]

    # Turn input filenames into a set of header objects
    if isinstance(inputs[0], str):
        hdrlist = [[],[],[],[]]
        for fname in inputs:
            #print 'Getting single template for : ',fname
            hdrs = getSingleTemplate(fname, extlist=extlist)
            for i in range(len(hdrs)): 
                hdrlist[i].append(hdrs[i])
    else:
        hdrlist = inputs
    num_files = len(inputs)
    
    # create list of combined PRIMARY/SCI headers for use in creating
    # the new table extensions
    tabhdrs = []
    for phdr,scihdr in zip(hdrlist[0],hdrlist[1]):
        tabhdrs.append(cat_headers(phdr,scihdr))

    # Determine what blending rules need to be merged to create the final 
    # blended headers. There will be a separate set of rules for each 
    # instrument, and all rules get merged into a composite set of rules that
    # get applied to all input headers regardless of instrument.
    # 
    # Instrument identification will be extracted from the INSTRUME keyword from
    # the PRIMARY header of each input
    #
    icache = {}
    for i in range(num_files):
        ph = hdrlist[0][i]
        inst = ph['instrume'].lower()
        hlist = [hdrlist[0][i],hdrlist[1][i],hdrlist[2][i],hdrlist[3][i]]
        if inst not in icache:
            # initialize the appropriate class for this data's instrument
            inst_class = KeywordRules(inst.lower())
            # Interpret rules for this class based on image that 
            # initialized this instrument's rules
            for h in hlist:
                inst_class.interpret_rules(h)

            # Now add this interpreted class to the cache
            icache[inst] = inst_class

    # Create final merged set of rules
    final_rules = None
    for inst in icache:
        if final_rules is None:
            final_rules = icache[inst]
        else:
            final_rules.merge(icache[inst])

    # Apply rules to each set of input headers 
    new_headers = []
    i=0
    for hdrs in hdrlist:
        newhdr, newtab = final_rules.apply(hdrs)
        new_headers.append(newhdr)
    # Create extension table from list of all combined PRI/SCI headers
    tabhdr, newtab = final_rules.apply(tabhdrs)

    if len(newtab) > 0:
        # Now merge the results for all the tables into a single table extension
        new_table = pyfits.new_table(newtab)
    else:
        new_table = None
    return new_headers,new_table
        
#### Classes for managing keyword rules
class KeywordRules(object):
    
    rules_name_suffix = '_header.rules'
    
    def __init__(self,instrument):
        """ Read in the rules used to interpret the keywords from the specified 
            instrument image header.
        """
        self.instrument = instrument
        self.new_header = None
        self.rules_file = self.get_filename()
        if self.rules_file is None:
            errmsg = textutil.textbox('ERROR:\n'+
                '    No valid rules file found for instrument: %s\n'%
                    (self.instrument))
            print errmsg
            raise ValueError

        rfile = open(self.rules_file)
        self.rule_specs = rfile.readlines()
        rfile.close()
        
        self.rules = []
        self.section_names = []
        self.delete_kws = []
        
    def get_filename(self):
        """ Return name of rules file to be used
        It will use a local copy if present, and use the installed version
        by default. Any local copy will take precendence over the default rules.
        """
        rules_name = self.instrument.lower()+self.rules_name_suffix
        if os.path.exists(rules_name):
            rules_file = rules_name
        else:
            rules_file = os.path.join(os.path.dirname(__file__),rules_name)
            if not os.path.exists(rules_file):
                rules_file = None
        return rules_file
    
    def interpret_rules(self,hdrs):
        """ Convert specifications for rules from rules file
            into specific rules for this header(instrument/detector)
            
            This allows for expansion rules to be applied to rules 
            from the rules files (such as any wildcards or section titles).
        """        
        if isinstance(hdrs, tuple):
            hdrs = list(hdrs)
        if not isinstance(hdrs, list):
            hdrs = [hdrs]

        # remove boolean keywords from headers as they will not 
        # produce valid values for inclusion in the pyfits.BinTableHDU
        #
        for hdr in hdrs:
            kw_del_list = []
            for kw in hdr.ascard:
                if isinstance(kw.value,np.bool) and (
                    kw.key not in fits_required_bool_kws):
                    kw_del_list.append(kw.key)
            for kw in kw_del_list[-1::-1]:
                del hdr[kw]

        # apply rules to 'cleaned' headers
        for hdr in hdrs:
            new_rules = []
            for rule in self.rule_specs:
                irule,section_name,delete_kws = interpret_line(rule,hdr)
                if irule is not None and len(irule) > 0:
                    for i in irule: # prevent multiple rules for a kw
                        if i not in new_rules: new_rules.append(i)
                # Keep track of section names, so they can be deleted from 
                # the final output header
                if section_name is not None:
                    self.section_names.append(section_name)
                if len(delete_kws) > 0:
                    self.delete_kws.extend(delete_kws)

            if len(self.rules) > 0:
                # merge new set of rules with previously interpreted rules
                self.merge(new_rules)
            else:
                self.rules = new_rules
                        
        
    def merge(self,kwrules):
        """
        Merge a new set of interpreted rules into the current set
        The new rules, kwrules, can either be a new class or a whole new
        set of rules (like those obtained from using self.interpret_rules with
        a new header).
        """
        if isinstance(kwrules, KeywordRules):
            kwrules = kwrules.rules
        # Determine what rules are specified in kwrules that 
        #    are NOT in self.rules
        delta = list(set(kwrules) - set(self.rules))

        # extend self.rules with additional rules
        self.rules.extend(delta)

    def apply(self,headers,tabhdu=False):
        """ For a full list of headers, apply the specified rules to 
            generate a dictionary of new values and a table using
            fitsblender. This method can be called separately for
            each type of header that needs to be generated; such as,
            Primary, SCI, WHT, CTX,...
                       
            This method returns the new header and summary table
            as pyfits.Header and numpy.ndarray masked array or 
            pyfits.binTableHDU objects
        """
        # Apply rules to headers
        fbdict,fbtab = blender.fitsblender(headers,self.rules)

        # Determine which keywords are included in the table but not 
        # the new dict(header). These will be removed from the output
        # header altogether
        tabcols = fbtab.dtype.names
        hdrkws = fbdict.keys()
        del_kws = list(set(tabcols) - set(hdrkws))

        # Start with a copy of the template as the new header 
        # This will define what keywords need to be updated, as the rules
        # and input headers often include headers for multiple extensions in
        # order to build the complete table for all the keywords in the file
        # in one run
        new_header = headers[0].copy()

        # Delete all keywords from copy that are being moved into the table
        for kw in del_kws:
            if kw in new_header: del new_header[kw]
        
        # Remove section names from output header(s)
        for name in self.section_names:
            for indx,kw in zip(range(len(new_header),0,-1),new_header.ascard[-1::-1]):
                if name in str(kw.value): 
                    del new_header[indx-1]
                continue
            
        # Delete keywords marked in rules file 
        for kw in self.delete_kws:
            if kw in new_header: del new_header[kw]
            
        # Apply updated/blended values into new header, but only those
        # keywords which are already present in the 'template' new header
        # this allows the rules to be used on all extensions at once yet 
        # update each extension separately without making copies of kws from
        # one extension to another.
        for kw in fbdict:
            if kw in new_header:
                new_header.update(kw,fbdict[kw],savecomment=True)

        # Create summary table 
        if len(tabcols) > 0:
            if tabhdu:
                new_table = pyfits.new_table(fbtab)
            else:
                new_table = fbtab
        else:
            new_table = None
        return new_header,new_table
            
    def index_of(self,kw):
        """ Reports the index of the specified kw
        """
        indx = []
        for r,i in zip(self.rules, range(len(self.rules))):
            if r[0] == kw: indx.append(i)
        return indx
    
#### Utility functions 
def interpret_line(line, hdr):
    """ Generate the rule(s) specified by the input line from the rules file 
    """
    # Initialize output values
    rules = []
    section_name = None
    delete_kws = []
    # Ignore comment lines in rules file
    if line[0] == '#' or len(line.strip()) == 0:
        return rules,section_name,delete_kws
    # clean up input lines
    line = line.strip('\n')
    
    # strip off any comment from the line before parsing the line
    if '#' in line: line = line[:line.rfind('#')].strip()
    
    # Parse the line
    use_kws = True
    if delete_command in line:
        use_kws = False
        line = line.replace(delete_command,'').lstrip()
        
    if '/' in line:
        section_name = line.split('/')[1].strip()
        kws = find_keywords_in_section(hdr, section_name)
        if kws is not None:
            for kw in kws:
                if use_kws:
                    rules.append((kw,kw))
                else:
                    delete_kws.append(kw)
    else:
        kwnames = line.split()
        if '*' in line:
            kws = hdr[kwnames[0]].keys()
            kws2 = kws
        else:
            kws = [kwnames[0]]
            if len(kwnames) > 1: 
                kws2 = [kwnames[1]]
            else:
                kws2 = kws
                
        lrule = None
        if len(kwnames) > 2:
            indx = line.index(kwnames[2])
            lrule = line[indx:].strip() # rule used for new value in header

        # Interpret short-hand rules using dict
        if lrule is not None and len(lrule) > 0:
            if lrule in blender_funcs:                        
                lrule = blender_funcs[lrule]
            else:
                lrule = None
            # build separate rule for each kw
            for kw,kw2 in zip(kws,kws2):
                if use_kws:
                    new_rule = (kw,kw2,lrule,"ignore")
                    if new_rule not in rules: rules.append(new_rule)
                else:
                    delete_kws.append(kw)
        else:                
            for kw,kw2 in zip(kws,kws2):
                if use_kws:
                    new_rule = (kw,kw2)
                    if new_rule not in rules: rules.append(new_rule)
                else:
                    delete_kws.append(kw)
                    
    return rules,section_name,delete_kws

def find_keywords_in_section(hdr,title):
    """ Return a list of keyword names from hdr identified in the section 
        with the specified section title.
    """
    # Indentify card indices of start and end of specified section
    sect_start = None
    sect_end = None
    for i,kw in enumerate(hdr.ascard):
        if sect_start is None:
            if title in str(hdr[i]):
                sect_start = i
        else:
            if '/' in str(hdr[i]) and hdr[i] not in ['N/A',' ','']:
                sect_end = i
                break
    if sect_end is None: sect_end = len(hdr)
    if sect_start is None:
        return None

    # Now, extract the keyword names from this section
    section_keys = hdr.ascard[sect_start+1:sect_end-1].keys()
    # remove any blank keywords
    while section_keys.count('') > 0:
        section_keys.remove('')
    
    return section_keys

def getSingleTemplate(fname, extlist=['SCI', 'ERR', 'DQ']):
    """
    # Obtain default headers for output file based on a single input file
    # (Copied from outputimage module.)
    #
    # NOTE: Returns 'pyfits.Header' objects, not HDU objects!
    #
    """
        
    if fname is None:
        raise ValueError('No data files for creating FITS output.')
    
    froot,fextn = fileutil.parseFilename(fname)
    if fextn is not None:
        fnum = fileutil.parseExtn(fextn)[1]
    ftemplate = fileutil.openImage(froot,mode='readonly')
    prihdr = ftemplate['PRIMARY'].header.copy()
    del prihdr['pcount']
    del prihdr['gcount']

    if fname.find('.fits') > 0 and len(ftemplate) > 1:

        # Setup which keyword we will use to select each
        # extension...
        _extkey = 'EXTNAME'

        defnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[0])
        #
        # Now, extract the headers necessary for output (as copies)
        # 1. Find the SCI extension in the template image
        # 2. Make a COPY of the extension header for use in new output file
        if fextn is None:
            extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[0])
        else:
            extnum = (extlist[0],fnum)
        scihdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())
        #scihdr.update('extver',1)
        extnum_sci = extnum

        # Extract the header for additional extensions
        if len(extlist) > 1 and extlist[1] not in [None,'',' ','INDEF','None']:
            if fextn is None:
                extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[1])
            else:
                # there may or may not be a second type of extension in the template
                count = 0
                for f in ftemplate:
                    if f.header.has_key('extname') and f.header['extname'] == extlist[1]:
                        count += 1
                if count > 0:
                    extnum = (extlist[1],fnum)
                else:
                    # Use science header for remaining headers
                    extnum = (extlist[0],fnum)
        else:
            extnum = extnum_sci

        errhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())
        #errhdr.update('extver',1)
        errhdr.update('bunit','UNITLESS')
        
        if len(extlist) > 2 and extlist[2] not in [None,'',' ','INDEF','None']:
        
            if fextn is None:
                extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[2])
            else:
                count = 0
                for f in ftemplate:
                    if f.header.has_key('extname') and f.header['extname'] == extlist[2]:
                        count += 1
                if count > 0:
                    extnum = (extlist[2],fnum)
                else:
                    # Use science header for remaining headers
                    extnum = (extlist[0],fnum)
        else:
            extnum = extnum_sci
            
        dqhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())
        #dqhdr.update('extver',1)
        dqhdr.update('bunit','UNITLESS')

    else:
        # Create default headers from scratch
        scihdr = None
        errhdr = None
        dqhdr = None

    ftemplate.close()
    del ftemplate

    return prihdr,scihdr,errhdr,dqhdr  

def merge_tables_by_cols(tables):
    """ 
    Merge all input tables provided as a list of np.ndarray objects into a
    single np.ndarray object 
    """
    # build new list of dtypes for all cols from all tables
    # However, skip duplicate columns from subsequent tables
    # Only the values from the last table will be kept for duplicate columns
    new_dtypes = [] 
    new_cnames = []
    num_rows = len(tables[0])
    print '#'*40
    print ''
    print 'New merge...'
    print ''
    print '#'*40
    for t in tables:
        for cname in t.dtype.names:
            if cname not in new_cnames:
                new_cnames.append(cname)
                new_dtypes.append((cname,t[cname].dtype.str))
        print '#'*40
        print t.dtype
        print '#'*40
    # create an empty table with the combined dtypes from all input tables
    new_table = np.zeros((num_rows,),dtype=new_dtypes)
    print new_dtypes
    # copy data column-by-column from each input to new output table 
    for t in tables:
        for cname in t.dtype.names:
            print 'CNAME: ',cname,' with dtype: ',t[cname].dtype,new_table[cname].dtype
            new_table[cname] = t[cname]

    return new_table

def remove_distortion_keywords(hdr):
    """
    Remove WCS distortion related keywords from the input header
    """
    distortion_kws = ['TDDALPHA','TDDBETA',
                        'D2IMEXT','D2IMERR',
                        'DGEOEXT','NPOLEXT']
                        
    # Remove any reference to TDD correction from
    #    distortion-corrected products
    # We also need to remove the D2IM* keywords so that HSTWCS/PyWCS
    # does not try to look for non-existent extensions
    for kw in distortion_kws:
        if kw in hdr:
            del hdr[kw]

    # Remove '-SIP' from CTYPE for output product
    if 'ctype1' in hdr and hdr['ctype1'].find('SIP') > -1:
            hdr.update('ctype1', hdr['ctype1'][:-4])
            hdr.update('ctype2',hdr['ctype2'][:-4])

    # Remove SIP coefficients from DRZ product
    for k in hdr.items():
        if (k[0][:2] in ['A_','B_']) or (k[0][:3] in ['IDC','SCD'] and k[0] != 'IDCTAB') or \
        (k[0][:6] in ['SCTYPE','SCRVAL','SNAXIS','SCRPIX']):
            del hdr[k[0]]
    # Remove paper IV related keywords related to the
    #   DGEO correction here
    for k in hdr.items():
        if (k[0][:2] == 'DP'):
            del hdr[k[0]+'*']
            del hdr[k[0]+'.*']
            del hdr[k[0]+'.*.*']
        if (k[0][:2] == 'CP'):
            del hdr[k[0]]
    
def cat_headers(hdr1,hdr2):
    """
    Create new pyfits.Header object from concatenating 2 input Headers
    """
    nhdr = hdr1.copy().ascardlist()
    nhdr += hdr2.ascardlist()
            
    return pyfits.Header(nhdr)

def extract_filenames_from_drz(drzfile):
    """ 
    Returns the list of filenames with extensions of input chips that was
    used to generate the drizzle product.
    """
    phdr = pyfits.getheader(drzfile)
    fnames = phdr['D0*DATA'].values()
    return fnames

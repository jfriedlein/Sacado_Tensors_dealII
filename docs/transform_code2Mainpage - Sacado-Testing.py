#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 18:50:30 2019

@author: johannes
"""
# useful commands to clear the workspace and clear the console (MATLAB: clear,clc) from https://stackoverflow.com/questions/54943710/code-to-clear-console-and-variables-in-spyder
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

# Transform normal code into Doxygen mainpage
dirMainpage = "/home/johannes/DEALII/MA-Code/Sacado-Testing/docs/mainpage-framework.h"
dirCode = "/home/johannes/DEALII/MA-Code/Sacado-Testing/Sacado_example.cc"
mergedFile = "/home/johannes/DEALII/MA-Code/Sacado-Testing/mainpage.h"

do_replace = False
with open(mergedFile, "w") as outputFile:
    with open(dirMainpage, "r") as mainpageDoc:
        for line in mainpageDoc:
            if ( "\code" in line ):
                #outputFile.write(line) # output this found "\code" from the mainpage file
                # modify the code an write it into the mainpage file
                code = True
                firstLine = True
                #emptyLine = 0
                with open(dirCode, "r") as codeDoc:
                    for lineCode in codeDoc:
                        first_nonspace_index = len(lineCode) - len(lineCode.lstrip())
                        first_N_nonspace_letters = lineCode[first_nonspace_index:first_nonspace_index+3]
                        codeLast = code
                        
                        # instead of /// we indicate a completely commented out line by use no space after //, e.g. "//std::cout ..." is commented out but not a comment like "// This is a comment"
                        #if ( first_N_nonspace_letters == "///" ): # this indicates some completely outcommented code, so we keep the comment and don't create a text from this
                        #    lineCode = lineCode[1:len(lineCode)] # on
                        #    code = True
                        if ( first_N_nonspace_letters == "// " ):
                            if ( len(lineCode)==3 ): # this would be an empty line marked as comment
                                lineCode = "\n"
                            else:
                                lineCode = lineCode[(first_nonspace_index+3):len(lineCode)]
                            code = False
                        elif ( first_N_nonspace_letters == "//\n" ): # when it is an entirely empty comment line (without blank after //)
                            lineCode = "\n"
                            code = False
                        elif ( first_N_nonspace_letters == "" ): # line is an empty line
                            lineCode =" \n"
                            code = codeLast     # save the type of the last line
                        else:
                            # lineCode is used unchanged
                            code = True
                         
                        if ( firstLine == True ):
                            if ( code == True ): # if the first line is a code then we have to print the "\code"-command
                                codeLast = False # to trigger the corresponding "elif" below
                            elif ( code == False ):
                                codeLast = False # so we will just print the line as it is
                            firstLine = False
                                                
                        #if ( emptyLine == False ):
                        if ( code==codeLast ):
                            outputFile.write(lineCode)
                        elif ( codeLast==False and code==True ):
                            outputFile.write("\code\n")
                            outputFile.write(lineCode)
                        elif ( codeLast == True and code == False ):
                            outputFile.write("\endcode\n")
                            outputFile.write(lineCode)
                            
                    if ( code == True ): # ensure that the blocks are closed at the end
                        outputFile.write("\endcode\n")
                do_replace = True
        
            if ( do_replace == False ):
                outputFile.write(line)
        
            if ( do_replace == True and ("\endcode" in line) ): # when you find the "\endcode"-command in the source mainpage then stop replacing things, but only for the next line because we don't want to add this \endcode-command to the outputFile
                do_replace = False
                

            

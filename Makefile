# CXFORM multi-purpose makefile
# See INSTALL.TXT for more information
#
# 2000/06/14: Ed Santiago - Original version
# 2003/09/12: Ryan Boller - Modified for multiple platforms/languages
# 2004/03/18: Ryan Boller - Modified for OS auto-detect, Mac compatibility
# 2009/11/25: Ryan Boller - Added compiler flags for 64-bit Mac platform

SRC = 	cxform.c
ASM = 	cxform.a
DEF = 	cxform.def


# NOTE: If using MS Windows, this Makefile is only useful for gcc/cygwin builds.
#    Also note that these builds will not be compatible with IDL 6.0+, as
#    they are based on Microsoft Visual C++ 7.0.  If you need to build for
#    IDL 6.0+, use the "make_CXFORM_MSVC.bat" file.  You will also need
#    to have either MSVC++ 7.0 (not free) or MSVC++ Toolkit 2003 (free)
#    installed.


# NOTE: User needs to modify appropriate directory for system-specific location of IDL
IDL_DIR_UNIX=/usr/local/rsi/idl
IDL_DIR_MAC=/Applications/exelis/idl
IDL_DIR_WIN=c:\dev\RSI\IDL63


IDL_EXT_UNIX=$(IDL_DIR_UNIX)/external
IDL_EXT_MAC=$(IDL_DIR_MAC)/external
IDL_EXT_WIN=$(IDL_DIR_WIN)\external
IDL_LIB_WIN=$(IDL_DIR_WIN)\bin\bin.x86


# Detect which OS is running and select appropriate flags
all:		dll
so:
	@echo "OS type detected: "`uname`
	@case `uname` in \
		"SunOS")  make sun-cxform.so \
				"CFLAGS=-fPIC -shared -I$(IDL_EXT_UNIX)" ;;\
		"Darwin") make mac-cxform.so \
				"CFLAGS=-fPIC -arch x86_64 -I$(IDL_EXT_MAC)" ;;\
		"Linux")  make linux-cxform.so \
				"CFLAGS=-fPIC -shared -I$(IDL_EXT_UNIX)" ;;\
		*) echo "This operating system is not supported -- use make dll if under MS Windows" ;;\
	esac

so-c:
	@echo "OS type detected: "`uname`
	@case `uname` in \
		"SunOS")  make sun-cxform-c.so \
				"CFLAGS=-fPIC -shared" ;;\
		"Darwin") make mac-cxform-c.so \
				"CFLAGS=-fPIC -arch x86_64" ;;\
		"Linux")  make linux-cxform-c.so \
				"CFLAGS=-fPIC -shared" ;;\
		*) echo "This operating system is not supported -- use make dll-c if under MS Windows" ;;\
	esac


dll:		
	make cxform.dll "CFLAGS=-Wall -shared -g -I$(IDL_EXT_WIN)"
dll-c:		
	make cxform-c.dll "CFLAGS=-Wall -shared -g"


# MAC-SPECIFIC
mac-cxform.so:	cxform-auto.o  cxform-manual.o  cxform-dlm.o 
	gcc -arch x86_64 -flat_namespace -undefined suppress -bundle -o cxform.so $^ 

mac-cxform-c.so: cxform-auto.o  cxform-manual.o
	gcc -arch x86_64 -flat_namespace -undefined suppress -bundle -o cxform-c.so $^


# LINUX-SPECIFIC
linux-cxform.so:  cxform-auto.o  cxform-manual.o  cxform-dlm.o 
	ld -G -lm -o cxform.so $^

linux-cxform-c.so: cxform-auto.o  cxform-manual.o
	ld -G -lm -o cxform-c.so $^


# SUN-SPECIFIC
sun-cxform.so:  cxform-auto.o  cxform-manual.o  cxform-dlm.o 
	ld -G -lm -o cxform.so $^

sun-cxform-c.so: cxform-auto.o  cxform-manual.o
	ld -G -lm -o cxform-c.so $^




# former UNIX build command
#cxform-c.so:	cxform-auto.o  cxform-manual.o
#	ld -G -lm -o $@ $^


# WINDOWS-SPECIFIC
cxform.dll:	cxform-auto.o  cxform-manual.o  cxform-dlm.o
	gcc -shared -W1,--enable-auto-image-base -o $@ -W1,--out-implib=$@.a $(DEF) $^ $(IDL_LIB_WIN)\idl.lib
	dlltool --def $(DEF) --dllname $@ --output-lib $(ASM)

cxform-c.dll: cxform-auto.o  cxform-manual.o
	gcc -shared $^ -o $@



#####################
# C testing routines
main.exe: cxform-c.dll main.o
	gcc $^ -o $@

tester.o:
	gcc -c -o tester.o tester.c

tester.exe: cxform-c.dll tester.o
	gcc $^ -o $@

interpolation-issues: interpolation-issues.c
	gcc interpolation-issues.c cxform-auto.c cxform-manual.c -o interpolation-issues
	./interpolation-issues 

# Create test results using revision b00f20 in main branch ('added-interfaces')
# In previous revisions, tester would not run ("Illegal instruction: 4") due
# to a string allocation error. Error occured on system with 'gcc --version' returning
#   Configured with: --prefix=/Library/Developer/CommandLineTools/usr --with-gxx-include-dir=/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/usr/include/c++/4.2.1
#   Apple clang version 12.0.0 (clang-1200.0.32.29)
#   Target: x86_64-apple-darwin19.5.0
test-ref: tester.c main.c
	git checkout b00f20 >> test-log/test-ref-b00f20.log 2>&1
	@gcc tester.c cxform-auto.c cxform-manual.c -o tester >> test-log/test-ref-b00f20.log 2>&1
	@./tester >> test-log/test-ref-b00f20.log 2>&1
	@mv test_results.txt test-log/test_results-b00f20.txt
	@make --always-make main >> test-log/test-ref-b00f20.log 2>&1
	@./main > test-log/main_results-b00f20.txt
	@git checkout added-interfaces >> test-log/test-ref-b00f20.log 2>&1

# Get current repository revision hash, which has a correction that caused
# IGRF coefficients after 2015 to be interpolated with the 2015 value and zero
# due to a commit by Jillian.
REV=$(shell git rev-parse --short HEAD)
# Create test results using current revision of main branch ('added-interfaces')
test: tester.c main.c
	@git checkout added-interfaces >> test-log/test-$(REV).log 2>&1
	@gcc tester.c cxform-auto.c cxform-manual.c -o tester >> test-log/test-$(REV).log 2>&1
	@./tester >> test-log/test-ref-$(REV).log 2>&1

	@mv test_results.txt test-log/test_results-$(REV).txt
	@diff test-log/test_results-$(REV).txt test-log/test_results-$(REV).txt && echo "PASS: test_results file matches reference version."

	@make --always-make main >> test-log/test-$(REV).log 2>&1
	@./main > test-log/main_results-$(REV).txt
	@diff test-log/main_results-$(REV).txt test-log/main_results-$(REV).txt && echo "PASS: main_results file matches reference version."


main: main.c cxform-auto.c cxform-manual.c
	@echo "OS type detected: "`uname`
	@case `uname` in \
		"SunOS")  make sun-main \
				"CFLAGS=-fPIC -shared" ;;\
		"Darwin") make mac-main \
				"CFLAGS=-fPIC" ;;\
		"Linux")  make linux-main \
				"CFLAGS=-fPIC -shared" ;;\
		*) echo "This operating system is not supported -- use make main.exe if under MS Windows" ;;\
	esac


sun-main: sun-cxform-c.so main.o
	gcc ./cxform-c.so main.o -lm -o main
#	gcc $^ -lm -o $@


mac-main:
	@echo "NOTE TO MAC USERS: There is currently a problem linking a C routine to the cxform-c.so shared library.  The alternative is to compile cxform directly with the routine, as done here"
	gcc main.c cxform-auto.c cxform-manual.c -o main

linux-main: linux-cxform-c.so main.o
	gcc ./cxform-c.so main.o -lm -o main



#######################
clean:
	rm -f *.so *.sl *.o *.obj *.a core a.out cxform.dll cxform-c.dll cxform.lib cxform.exp main tester interpolation-issues main.exe tester.exe


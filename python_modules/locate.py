import socket
import getpass
def locate():
	err=0
	if "jasmin" in (socket.gethostname()):
		loc='/group_workspaces/jasmin2/ukca/vol2/'+getpass.getuser()+'/'
	elif "jc.rl" in (socket.gethostname()):
		loc='/group_workspaces/jasmin2/ukca/vol2/'+getpass.getuser()+'/'
	elif "atm" in (socket.gethostname()):
		loc='/scratch/'+getpass.getuser()+'/netscratch/um/'
	else:
		err=err+1
	if err>0:
		print socket.gethostname()
		print "UNKNOWN ANALYSIS PLATFORM - PLEASE DEFINE IN LOCATE"	
	return loc
def locate2():
        err=0
        if "jasmin" in (socket.gethostname()):
                loc='jasmin'
        elif "jc.rl" in (socket.gethostname()):
                loc='jasmin'
        elif "atm" in (socket.gethostname()):
                loc='cambridge'
        else:
                err=err+1
        if err>0:
                print socket.gethostname()
                print "UNKNOWN ANALYSIS PLATFORM - PLEASE DEFINE IN LOCATE"
        return loc


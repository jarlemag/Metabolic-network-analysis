#Setlocalvalues.py

    
for key in bigdict:
	if (key in locals()) == False:
		locals()[key] = dig[key]

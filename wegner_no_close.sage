import wegner as wg



li = range(77)
del li[63]#cxa4x5555
del li[62]#cxa4x555x4
del li[53]#c8a3x55x55
del li[50]#c7a55x5
del li[22]#cxa555555






for i in li:
    s = wg.HRCs[i]
    wg.check_all_configuration_realizations_temp2(s)
    print ("-"*60+"\n")*5
print "Mission complete!"

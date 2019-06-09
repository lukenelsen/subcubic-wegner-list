
from checking_realizations import *

NewTargetSet = TargetSet[:]
for i in range(len(NewTargetSet)):
    rc = NewTargetSet[i]
    if len(rc) == 4:  # c3a3, c3a4, c3a6, c4a4
        NewTargetSet[i] = "cxa"+rc[1]+rc[3]

DataDict = {x:{} for x in NewTargetSet}

# First, add the names of the configurations in paper format.
DataDict['cxa33']['paper name'] = r'33'

DataDict['cxa3x3']['paper name'] = r'3\*3'

DataDict['cxa3xx3']['paper name'] = r'3\*\*3'

DataDict['cxa3xxx3']['paper name'] = r'3\*\*\*3'

DataDict['cxa3xx5x3']['paper name'] = r'3\*\*5\*3'

DataDict['cxa34']['paper name'] = r'34'

DataDict['cxa3x4']['paper name'] = r'3\*4'

DataDict['cxa35x4']['paper name'] = r'35\*4'

DataDict['cxa3x54']['paper name'] = r'3\*54'

DataDict['cxa3xx54']['paper name'] = r'3\*\*54'

DataDict['cxa3xxx54']['paper name'] = r'3\*\*\*54'

DataDict['cxa355']['paper name'] = r'355'

DataDict['cxa375']['paper name'] = r'375'

DataDict['cxa3x555']['paper name'] = r'3\*555'

DataDict['cxa36']['paper name'] = r'36'

DataDict['cxa356']['paper name'] = r'356'

DataDict['cxa44']['paper name'] = r'44'

DataDict['cxa454']['paper name'] = r'454'

DataDict['cxa464']['paper name'] = r'464'

DataDict['cxa474']['paper name'] = r'474'

DataDict['cxa4x54']['paper name'] = r'4\*54'

DataDict['cxa4x55x4']['paper name'] = r'4\*55\*4'

DataDict['cxa4x555x4']['paper name'] = r'4\*555\*4'

DataDict['cxa4x4x4x4']['paper name'] = r'4\*4\*4\*4'

DataDict['cxa4x45']['paper name'] = r'4\*45'

DataDict['cxa4555']['paper name'] = r'4555'

DataDict['cxa4x5555']['paper name'] = r'4\*5555'

DataDict['cxa456']['paper name'] = r'456'

DataDict['cxa535']['paper name'] = r'535'

DataDict['cxa555555']['paper name'] = r'555555'

DataDict['cxa546']['paper name'] = r'546'

DataDict['c3a57']['paper name'] = r'3:57'

DataDict['c4a55']['paper name'] = r'4:55'

DataDict['c4a56']['paper name'] = r'4:56'

DataDict['c4a57']['paper name'] = r'4:57'

DataDict['c4a66']['paper name'] = r'4:66'

DataDict['c4a585']['paper name'] = r'4:585'

DataDict['c4a676']['paper name'] = r'4:676'

DataDict['c4a686']['paper name'] = r'4:686'

DataDict['c5a555']['paper name'] = r'5:555'

DataDict['c5a556']['paper name'] = r'5:556'

DataDict['c5a565']['paper name'] = r'5:565'

DataDict['c5a566']['paper name'] = r'5:566'

DataDict['c5a575']['paper name'] = r'5:575'

DataDict['c5a656']['paper name'] = r'5:656'

DataDict['c5a666']['paper name'] = r'5:666'

DataDict['c5a4x55']['paper name'] = r'5:4\*55'

DataDict['c5a5x65']['paper name'] = r'5:5\*65'

DataDict['c5a5585']['paper name'] = r'5:5585'

DataDict['c5a5x66']['paper name'] = r'5:5\*66'

DataDict['c5a55x6']['paper name'] = r'5:55\*6'

DataDict['c5a56x6']['paper name'] = r'5:56\*6'

DataDict['c5a6x66']['paper name'] = r'5:6\*66'

DataDict['c7a3xx4']['paper name'] = r'7:3\*\*4'

DataDict['c7a3xx5']['paper name'] = r'7:3\*\*5'

DataDict['c7a4xx4']['paper name'] = r'7:4\*\*4'

DataDict['c7a4x55']['paper name'] = r'7:4\*55'

DataDict['c7a4xx55']['paper name'] = r'7:4\*\*55'

DataDict['c7a4x5xx5']['paper name'] = r'7:4\*5\*\*5'

DataDict['c7a55x5']['paper name'] = r'7:55\*5'

DataDict['c8a3xx4']['paper name'] = r'8:3\*\*4'

DataDict['c8a3xxx4']['paper name'] = r'8:3\*\*\*4'

DataDict['c8a35x5']['paper name'] = r'8:35\*5'

DataDict['c8a35xx5']['paper name'] = r'8:35\*\*5'

DataDict['c8a35xxx5']['paper name'] = r'8:35\*\*\*5'

DataDict['c8a35xxxx5']['paper name'] = r'8:35\*\*\*\*5'

DataDict['c8a3x55x55']['paper name'] = r'8:3\*55\*55'

DataDict['c8a45xx4']['paper name'] = r'8:45\*\*4'

DataDict['c8a455x5']['paper name'] = r'8:455\*5'

DataDict['c8a4x5x4x4']['paper name'] = r'8:4\*5\*4\*4'

DataDict['c9a3xx4']['paper name'] = r'9:3\*\*4'

DataDict['c9a3xxx4']['paper name'] = r'9:3\*\*\*4'

DataDict['c9a35x5']['paper name'] = r'9:35\*5'

DataDict['c9a455x4']['paper name'] = r'9:455\*4'

DataDict['c9a455x5']['paper name'] = r'9:455\*5'

DataDict['c9a545x5']['paper name'] = r'9:545\*5'



# Next, pull relevant data from file:  numbers for partitions, core subgraphs, and realizations, and also runtime.
DataDict['total'] = {s:0 for s in ['partition count','core subgraph count','realization count','runtime']}

# Partitions
pfile = open('partitions.txt','r')
for line in pfile:
    #First, get the file index since they are out of order.  Depending on our preprocessing, we might have RC_##: or RC_##.out, and either RC_0# or RC_# for single digits.
    #They probably start at 1 and go up to 76, but this should be verified.
    if line[4] in ['.',':']:  # One version of single digit.  The other version is captured with the double-digit cases.
        i = line[3]
    else:
        i = line[3:5]
    num = line[line.index('    ')+4:-1]  # The -1 is to cut off the new line character.  #partitions.txt is four spaces, the other three are three spaces.
    DataDict['total']['partition count'] += int(num)
    # Now add commas.
    q = (len(num)-1)/3  # Should be int.
    r = len(num)%3
    if r == 0:
        r = 3
    pnum = num[:r]
    a = 0
    while a < q:
        pnum += ','+num[r+3*a:r+3*a+3]
        a += 1
    DataDict[NewTargetSet[int(i)-1]]['partition count'] = pnum
pfile.close()

# Core Subgraphs
pfile = open('core_subgraphs.txt','r')
for line in pfile:
    #First, get the file index since they are out of order.  Depending on our preprocessing, we might have RC_##: or RC_##.out, and either RC_0# or RC_# for single digits.
    #They probably start at 1 and go up to 76, but this should be verified.
    if line[4] in ['.',':']:  # One version of single digit.  The other version is captured with the double-digit cases.
        i = line[3]
    else:
        i = line[3:5]
    num = line[line.index('   ')+3:-1]  # The -1 is to cut off the new line character.  #partitions.txt is four spaces, the other three are three spaces.
    DataDict['total']['core subgraph count'] += int(num)
    # Now add commas.
    q = (len(num)-1)/3  # Should be int.
    r = len(num)%3
    if r == 0:
        r = 3
    pnum = num[:r]
    a = 0
    while a < q:
        pnum += ','+num[r+3*a:r+3*a+3]
        a += 1
    DataDict[NewTargetSet[int(i)-1]]['core subgraph count'] = pnum
pfile.close()
pfile.close()

# Realizations
pfile = open('realizations.txt','r')
for line in pfile:
    #First, get the file index since they are out of order.  Depending on our preprocessing, we might have RC_##: or RC_##.out, and either RC_0# or RC_# for single digits.
    #They probably start at 1 and go up to 76, but this should be verified.
    if line[4] in ['.',':']:  # One version of single digit.  The other version is captured with the double-digit cases.
        i = line[3]
    else:
        i = line[3:5]
    num = line[line.index('   ')+3:-1]  # The -1 is to cut off the new line character.  #partitions.txt is four spaces, the other three are three spaces.
    DataDict['total']['realization count'] += int(num)
    # Now add commas.
    q = (len(num)-1)/3  # Should be int.
    r = len(num)%3
    if r == 0:
        r = 3
    pnum = num[:r]
    a = 0
    while a < q:
        pnum += ','+num[r+3*a:r+3*a+3]
        a += 1
    DataDict[NewTargetSet[int(i)-1]]['realization count'] = pnum
pfile.close()


def timestring_int(total):
    # total is a float counting seconds
    # Displays time length:  #d #h #m #s, only showing the ones more than 0.  If s < 1, then displays "< 1s".
    
    t = int(total)  # cut off the noninteger part
    print_str = ''
    d = t/86400
    t %= 86400
    if d:
        print_str += str(d)+'d '
    h = t/3600
    t %= 3600
    if h or d:
        print_str += str(h)+'h '
    m = t/60
    t %= 60
    if m or h or d:
        print_str += str(m)+'m '
    if t or m or h or d:
        print_str += str(t)+'s'
    else:
        print_str = r'$<$ 1s'
    return print_str


# Runtimes
pfile = open('runtimes.txt','r')
for line in pfile:
    #First, get the file index since they are out of order.  Depending on our preprocessing, we might have RC_##: or RC_##.out, and either RC_0# or RC_# for single digits.
    #They probably start at 1 and go up to 76, but this should be verified.
    if line[4] in ['.',':']:  # One version of single digit.  The other version is captured with the double-digit cases.
        i = line[3]
    else:
        i = line[3:5]
    m = line[line.index('   ')+3:line.index('m ')]  # minutes
    s = line[line.index('m ')+2:-2]  # seconds (-2 takes care of new line character and 's')
    sec = 60*float(int(m))+float(int(s[:-2]))+float(int(s[-1])/10.0)
    DataDict['total']['runtime'] += sec
    DataDict[NewTargetSet[int(i)-1]]['runtime'] = timestring_int(sec)
pfile.close()



#print DataDict['total']















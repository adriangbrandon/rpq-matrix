import pandas as pd


M = {
    "v /* c": 1,
    "v * c": 2,
    "v + c": 3,
    "c * v": 4,
    "c /* v": 5,
    "v / c": 6,
    "v / v": 7,
    "v */* c": 8,
    "v |* c": 9,
    "v | v": 10,
    "v ^ v": 11,
    "v /* v": 12,
    "v */*/*/*/* c": 13,
    "v * v": 14,
    "v /? c": 15,
    "v + v": 16,
    "v /+ c": 17,
    "v ^/ c": 18,
    "v || v": 19,
    "v /^ v": 20,
    "v | c": 21
}


def get_dataframe():
	queries = []
	# Open the file in read mode
	with open('paths.tsv', 'r') as file:
	    # Read each line using a loop
	    for line in file:
	        # Process the line (print it in this example)
	       queries.append(line.strip())

	types = []
	# Open the file in read mode
	with open('queries-classified.tsv', 'r') as file:
	    # Read each line using a loop
	    for line in file:
	        # Process the line (print it in this example)
	    	t = line.strip()
	    	if t in M:
	    		types.append(M[t])
	    	else:
	    		types.append(0)

	negated = []
	for q in queries:
		negated.append("%" in q)

	s= set()
	with open('queries-selected.tsv', 'r') as file:
	    # Read each line using a loop
	    for line in file:
	        s.add(int(line.strip()))

	selected = []
	for i in range(1, len(queries)+1):
		if i in s:
			selected.append(True)
		else:
			selected.append(False)

	data = {"Id": range(1, len(queries)+1), "Q": queries, "Type": types, "Negated" : negated, "Selected": selected}
	df = pd.DataFrame(data)
	return df

def get_selected_queries(t):
	return df[(df['Type'] == t) & (df['Selected'] == True)]

def get_queries_to_remove(data):
	to_remove=[]
	for index, row in data.iterrows():
		if row['Negated']:
			to_remove.add(index)
	return to_remove



df = get_dataframe()

print(df[df['Negated']==True])

for i in range(1, 22):
	a = get_selected_queries(1)
	rem = get_queries_to_remove(a)
	out = open('type' + str(i) + 'remove.tsv', 'w')
	for r in rem:
		our.write(str(r) + '\n')
	out.close()


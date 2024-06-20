import copy


array1 = [2,4,5,6 ]
array2 = copy.deepcopy(array1)

print(array1)
print(array2)

array2[1] = 1000

print(array1)


# mx = array1[1]

# array1[1] = 200
# print(mx)


import numpy as np

A = np.array([[1, 1, 0], [1, 0, 1]])
B = np.zeros([2,3]) + 5

B[A >= 1] = 2
#print(A)

#B = A[A == 1]
# print(B.shape)
a, b = B.shape
# print(a, b)
# def test_function():
#     init_test.A = np.append(init_test.A, 1)

RH = 40
temp = 60

# np.save(f'RH {RH}, temp {temp}', A)

# C = np.load(f'RH {RH}, temp {temp}.npy')

print('Do you want to save stuff Y/N?')
answer = input()
if answer == 'Y':
    print(answer)
# print(A)
# print(C)
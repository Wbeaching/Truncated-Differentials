

import numpy as np
import math
import time

#MC_matrix = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]

def MC(A):
      return [A[1]^A[2]^A[3], A[0]^A[2]^A[3], A[0]^A[1]^A[3], A[0]^A[1]^A[2]]


def ham(n):
      ans=0
      if n==0:
            return 0
      while n>0:
            n=n&(n-1)
            ans+=1
      return ans

def list_to_num(A):
      s = ''
      for i in range(0,len(A)):
            s = s + str(A[i])
      return int(s,2) 


def code_diff_MC_out(A):
      for i in range(0,len(A)):
            if A[i] > 0:
                  A[i] = 1
      return A


def fixed_input_column(input_num):
      bin_input = bin(input_num)[2:].zfill(4)
      list_input = []
      for i in range(0,4):
            list_input.append(int(bin_input[i]))
      print(list_input)

      MC_out = []
      
      for a in range(list_input[0], 15* list_input[0] + 1 ):
            for b in range(list_input[1],15 * list_input[1] + 1): 
                  for c in range(list_input[2], 15 * list_input[2] + 1):
                        for d in range(list_input[3],15 * list_input[3] + 1):
                              MC_out.append(code_diff_MC_out(MC([a,b,c,d])))

      #print(MC_out)
      #print([a,b,c,d])
      N_out = [0 for i in range(16)]

      for i in range(len(MC_out)):
            N_out[list_to_num(MC_out[i])] += 1

      #print(sum(N_out))
      Pr_out = []
      for i in range(16):
            Pr_out.append(N_out[i]/sum(N_out))

      return Pr_out
      

def column_diff_propagation_pro_table():
      T = []
      for i in range(16):
            T.append(fixed_input_column(i))
      return T

T_col = column_diff_propagation_pro_table()

def fixed_input_state(input_num  ):

      #T_col = column_diff_propagation_pro_table()
      
      int_str = bin(input_num)[2:].zfill(16)
      
      col_1_num = int(int_str[0] + int_str[10] + int_str[5] + int_str[15], 2)
      col_2_num = int(int_str[11] + int_str[1] + int_str[14] + int_str[4], 2)
      col_3_num = int(int_str[6] + int_str[12] + int_str[3] + int_str[9], 2)
      col_4_num = int(int_str[13] + int_str[7] + int_str[8] + int_str[2], 2)

      #print(T_col)
      
      pr_col_1 = T_col[col_1_num]
      pr_col_2 = T_col[col_2_num]
      pr_col_3 = T_col[col_3_num]
      pr_col_4 = T_col[col_4_num]
      
      Pr_state = [0 for j in range(0, 2**16)]

      for output_num in range(0,2**16):
            output_str = bin(output_num)[2:].zfill(16)
            out_col_1_num = int(output_str[0] + output_str[4] + output_str[8] + output_str[12], 2)
            out_col_2_num = int(output_str[1] + output_str[5] + output_str[9] + output_str[13], 2)
            out_col_3_num = int(output_str[2] + output_str[6] + output_str[10] + output_str[14], 2)
            out_col_4_num = int(output_str[3] + output_str[7] + output_str[11] + output_str[15], 2)

            Pr_state[output_num] = pr_col_1[out_col_1_num] * pr_col_2[out_col_2_num] *  pr_col_3[out_col_3_num] * pr_col_4[out_col_4_num] 
      
      return Pr_state


def cal_str_diff_pro(diff_in, diff_out, T_col = T_col):
      print(fixed_input_state(int(diff_in,2))[int(diff_out,2)])
      return math.log(fixed_input_state(int(diff_in,2))[int(diff_out,2)],2)


def Two_Round_Test(diff_in, diff_out, T_col = T_col):
      pr = 0
      diff_one_out = fixed_input_state(int(diff_in, 2))
      for i in range(0, 2**16):
            if diff_one_out[i] >= 2**-10:
                  pr = pr + diff_one_out[i] * fixed_input_state(i)[int(diff_out,2)]
            
      return math.log(pr,2)


def get_TPT():
      TPT = []
      for i in range(0,2**16):
            TPT.append(fixed_input_state(i))
            #print(fixed_input(i))
      return TPT


def Mul_TPT(log_r):
      
      start = time.time()
      A = get_TPT()
      end = time.time()
      print("contrust matrix: ",  end - start)

      start = time.time()
      for i in range(log_r):
            A = np.dot(A,A)
      end = time.time()
      print("matrix mul",  end - start)
      return A


#random permutation j-th column
def ran_pr(j , cells ):
      
      pr = (1/(1 - 2**(-64))) * (1/16)**(cells - ham(j)) * (15/16)**ham(j)
      if j == 0:
            pr = 0
      #print(pr)
      #print(math.log(pr,2))
      return pr

def ran_matrix_row(cells ):
      A = []
      for i in range(2**cells):
            A.append(ran_pr(i,cells))
      #np.savetxt('ran_matrix.txt', A, delimiter=' ') 
      return A


def Find_high_pro(T ,R, cells ):
      #Delat_T = np.array([0 for i in range(2**cells)])
      index = []
      #zero_index = []
      for i in range(1,2**cells):
            for j in range(1,2**cells):
                  #if (T[i][j] == 0):
                        #zero_index.append((i,j))
                  if T[i][j]  >= 2**1 * R[j]:  #or (2*ran_matrix_row(cells)[j]) / T[i][j] <= 2**(-0.1):
                        index.append((i,j, math.log(T[i,j],2),math.log(R[j],2)))

      return index


print("GO! ")
start = time.time()

cells = 16

print('T4:')
T = Mul_TPT(2)  #T^4
print('mul')
R = ran_matrix_row(cells)
np.savetxt('Find_high_pro_Midori_4r.txt', Find_high_pro(T ,R, cells ), delimiter=' ')

print('T5:')
T = np.dot(T, Mul_TPT(0))  #T^5
np.savetxt('Find_high_pro_Midori_5r.txt', Find_high_pro(T ,R, cells ), delimiter=' ')

print('T6:')
T = np.dot(T, Mul_TPT(0))  #T^6
np.savetxt('Find_high_pro_Midori_6r.txt', Find_high_pro(T ,R, cells ), delimiter=' ')

print('T7:')
T = np.dot(T, Mul_TPT(0))  #T^7
np.savetxt('Find_high_pro_Midori_7r.txt', Find_high_pro(T ,R, cells ), delimiter=' ')

print('T8:')
T = np.dot(T, Mul_TPT(0))  
np.savetxt('Find_high_pro_Midori_8r.txt', Find_high_pro(T ,R, cells ), delimiter=' ')

print("END ")
end = time.time()
print('time: ', end - start)







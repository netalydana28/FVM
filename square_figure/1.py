#Write a Python program to write a list to a file.

# example_list = ["", "2", "3", "4"]

list_ind = [19,18]

x0, x1 =52,60
row_range = 10

with open('input.txt', 'w') as f:
    # for m in range(7):
    for i in range(x0,x1):
        # for j in range(1,4):
        f.write(f'{i},{list_ind[0]},{list_ind[1]},{list_ind[2]},{list_ind[3]}'+'\n')
        
        # list_ind[0]+=1
        # list_ind[1]=1
        # x0+=19
        # x1+=19

        
        # if row_range == 10:
        #     row_range = 9
        #     x1+=row_range
        # elif row_range == 9:
        #     row_range = 10
        #     x1+=row_range
        
        list_ind[0]+=1
        list_ind[1]+=1
        list_ind[2]+=1
        list_ind[3]+=1


        

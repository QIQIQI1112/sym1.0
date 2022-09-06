from rigid import *
from symdata import *
from itertools import combinations, combinations_with_replacement

# works number == 2
# ratio 1:1
def divide_list_to_equal_sum(list_input, ratio):
    # print(f'list_input:{list_input}')
    # print(f'ratio:{ratio}')
    list_len = len(list_input)
    check = 0
    for i in range(1, list_len):
        comb = list(combinations(list_input, i))
        # print(f'comb:{comb}')
        comb_len = len(comb)
        for j in range(comb_len):
            list_a = list(comb[j])
            list_inter = list_a.copy()
            list_b = [n for n in list_input if not n in list_inter or list_inter.remove(n)]
            # print(list_a)
            # print(list_b)
            sum_a = 0
            sum_b = 0
            for k in range(len(list_a)):
                sum_a += float(list_a[k][1])
            for k in range(len(list_b)):
                sum_b += float(list_b[k][1])
            ratio_check = round(abs((float(sum_a / sum_b)) - float((ratio[1] / ratio[2]))), 4)
            # print(check_4)
            # print('-------------------')
            if ratio_check < 0.001:
                check = 1
            if check == 1:
                break
        if check == 1:
            break
    if check == 1:
        return 'yes'
    else:
        return 'no'


def divide_list_to_equal_sum_count(list_input, ratio):
    # print(f'list_input:{list_input}')
    # print(f'ratio:{ratio}')
    case = 0
    list_len = len(list_input)
    check = 0
    for i in range(1, list_len):
        comb = list(combinations(list_input, i))
        # print(f'comb:{comb}')
        comb_len = len(comb)
        for j in range(comb_len):
            list_a = list(comb[j])
            list_inter = list_a.copy()
            list_b = [n for n in list_input if not n in list_inter or list_inter.remove(n)]
            # print(list_a)
            # print(list_b)
            sum_a = 0
            sum_b = 0
            for k in range(len(list_a)):
                sum_a += float(list_a[k][1])
            for k in range(len(list_b)):
                sum_b += float(list_b[k][1])
            ratio_check = round(abs((float(sum_a / sum_b)) - float((ratio[1] / ratio[2]))), 4)
            # print(check_4)
            # print('-------------------')
            if ratio_check < 0.001:
                case += 1
    return case


def rigid_wyckoff_sym(rigid_type, sym_no):
        if rigid_type == 'Tetrahedron':
            wyckoff_sym = Tetrahedron(sym_no, 2.0).wyckoff_sym()
        return wyckoff_sym


class SymCases:
    """
    ---------------------------------------------------------------
    todo list:
    mutiple rigid body before sym transform
    mutiple type of rigid bodies
    combine case_count and case_list
    ---------------------------------------------------------------
    ratio = array
    rigid_at_general = 1: rigid body allow to be at general position
                     = 0: not allowd to be at general position
    rigid_number: the max number of rigid bodies
    eg: sym = SymCases(tetrahedron, np.array([1, 3]), 1, 3)
    """    
    def __init__(self, rigid_type, ratio, rigid_at_general, rigid_number):
        self.rigid_type = rigid_type
        self.ratio = ratio
        self.gen = rigid_at_general
        self.rigid_no = rigid_number

    """
    mode = 1: count all number no matter how large
         = 0: TBD if number is large than 100000
    """
    def case_count(self, no_low, no_high, mode):
        dic_count = {}
        # single atom type number
        atom_type_no = len(self.ratio) - 1
        for sym_no in range(no_low, no_high + 1):
            case_dic = {}
            case_no = 0
            sg_len = len(sg[f'sg_{sym_no}'])
            # generate muti_list
            muti_list = []
            for m in range(1, sg_len):
                muti_list.append([f's{m}', int(sg[f'sg_{sym_no}'][f's{m}'][0])])
                # muti_list.append(int(sg[f'sg_{sym_no}'][f's{m}'][0]))

            # rigid_at_general = 0
            case_count = 0
            case_count_final = 0
            sg_range = 0
            if self.gen == 0:
                sg_range = sg_len
            elif self.gen == 1:
                sg_range = sg_len + 1
            # print(f'sg_range:{sg_range}')
            wyckoff_sym = rigid_wyckoff_sym(self.rigid_type, sym_no)
            for i in range(1, sg_range):
                if sg[f'sg_{sym_no}'][f's{i}'][1] in wyckoff_sym:
                    site_muti = int(sg[f'sg_{sym_no}'][f's{i}'][0])
                    if site_muti <= self.rigid_no:
                        ratio_total = 0
                        for j in range(1, len(self.ratio)):
                            ratio_total += self.ratio[j]
                        single_atom_no = site_muti * ratio_total
                        # the max number of atoms can be at general position
                        gen_max = int(single_atom_no / int(sg[f'sg_{sym_no}'][f's{sg_len}'][0]))
                        gen_muti = sg[f'sg_{sym_no}'][f's{sg_len}'][0]

                        if sg[f'sg_{sym_no}'][f's{i}'][3] == 1:
                            print(f'muti_list:{muti_list}')
                            for gen_no in range(gen_max + 1):
                                # for n in range(0, len(muti_list) + 1):
                                for n in range(0, single_atom_no):
                                    count2 = 0
                                    enum = list(combinations_with_replacement(muti_list, n))
                                    print(f'enum:{enum}')
                                    enum_inter = enum.copy()
                                    # print(enum_inter)
                                    if mode == 0 and len(enum) > 100000:
                                        case_count = 'TBD'
                                        break
                                    elif mode == 1 or (mode == 0 and len(enum) <= 100000):
                                        if len(enum) != 0:
                                            len_k = len(enum[0])
                                            for k in enum:
                                                # print(f'k:{k}')
                                                sum_k = 0
                                                for k1 in range(len_k):
                                                    sum_k += int(k[k1][1])
                                                if sum_k + (gen_muti * gen_no) != single_atom_no:
                                                    enum_inter.remove(k)
                                        # print(f'rigid_position:{i}')
                                        # print(f'atom_at_general:{gen_no}')
                                        # print(f'n:{n}')
                                        # print(enum)
                                        # print(enum_inter)
                                        # print(f'sum:{sum_k + (gen_muti * gen_no)}, single_atom_no:{single_atom_no}')
                                        # print('-----------')
                                        if len(enum_inter) != 0:
                                            case_no += 1
                                            case_list = []
                                            for k2 in range(len(enum_inter)):
                                                case_list_inter = []
                                                for k3 in range(len_k):
                                                    case_list_inter.append(enum_inter[k2][k3][0])
                                                for k4 in range(gen_no):
                                                    case_list_inter.append(f's{sg_len}')
                                                case_list.append(case_list_inter)
                                            case_dic[f'{case_no}'] = [f'{i}', case_list]

                                            count2 += 1
                                            # print(f'sg:{sym_no}')
                                            # print(f'site:s{i}')
                                            # print(f'gen number:{gen_no}')
                                            # print(f'atom at special:{n}')
                                            # print(f'expected:{len(enum)}')
                                            # print(f'current:{count2}')
                                            # print('-------------------')
                                        enum_len = len(enum_inter)
                                        print(f'enum_inter:{enum_inter}')
                                        print(f'enmu_len:{enum_len}')
                                        if atom_type_no == 1:
                                            case_count += enum_len
                                            # print(f'case count:{case_count}')
                                        elif atom_type_no == 2:
                                            print(f'enmu_inter:{enum_inter}')
                                            for k in range(enum_len):
                                                enum_list = list(enum[k])
                                                print(f'enmu_list:{enum_list}')
                                                check_atom = divide_list_to_equal_sum(enum_list, self.ratio)
                                                if check_atom == 'yes':
                                                    case_count += 1
                                if case_count == 'TBD':
                                    break
                            if case_count == 'TBD':
                                break
            print(f'case_count:{case_count}')
            print(f'case_dic:{case_dic}')
            if atom_type_no == 1:
                pass
            elif atom_type_no == 2:
                for n in range(1, len(case_dic) + 1):
                    dic_check_list = case_dic[f'{n}'][1]
                    dic_remove_list = []
                    for m in range(len(case_dic[f'{n}'][1])):
                        wyckoff_list = [case_dic[f'{n}'][1][m]]
                        check_list = []
                        for l in range(len(case_dic[f'{n}'][1][m])):
                            wyckoff_symble = case_dic[f'{n}'][1][m][l]
                            wyckoff_muti = sg[f'sg_{sym_no}'][wyckoff_symble][0]
                            check_add = [wyckoff_symble, wyckoff_muti]
                            check_list.append(check_add)
                            # print(f'check_add:{check_add}')
                        print(f'check_list:{check_list}')
                        check_atom = divide_list_to_equal_sum(check_list, self.ratio)
                        if check_atom == 'yes':
                            pass
                        elif check_atom == 'no':
                            dic_remove_list.append(wyckoff_list)
                    print(f'remove_list:{dic_remove_list}')
                    for p in dic_remove_list:
                        dic_check_list.remove(p)
                    
            # get rid of atom at fixed wyckoff position
            '''
            (list become empty)
            '''
            s_list = []
            for i in range(sg_len):
                s_x = sg[f'sg_{sym_no}'][f's{i + 1}'][2][0]
                s_y = sg[f'sg_{sym_no}'][f's{i + 1}'][2][1]
                s_z = sg[f'sg_{sym_no}'][f's{i + 1}'][2][2]
                if s_x != 'x' and s_y != 'y' and s_z != 'z':
                    s_list.append(f's{i + 1}')

            s_list_len = len(s_list)
            # print(s_list)
            if s_list_len != 0:
                dic_len = len(case_dic)
                for i in range(dic_len):
                    check_list = case_dic[f'{i + 1}'][1]
                    # print(f'check_list:{check_list}')
                    # print(f'check_list_len:{len(check_list)}')
                    j_remove = []
                    for j in range(len(check_list)):
                        s_check = 0
                        for k in range(s_list_len):
                            if check_list[j].count(s_list[k]) >= 2:
                                # print(f'j:{j}')
                                # print(f'{s_list[k]} count:{check_list[j].count(s_list[k])}')
                                s_check = 1
                                # print(f's_check:{s_check}')
                                # print('------------------------------')
                                break
                        if s_check == 1:
                            j_remove.append(j)
                    j_remove_len = len(j_remove)
                    print(f'j_remove:{j_remove}')
                    j_remove_list = []
                    if j_remove_len != 0:
                        for n in j_remove:
                            j_remove_list.append(check_list[n])
                    for m in j_remove_list:
                        check_list.remove(m)

            # count number
            if atom_type_no == 1:
                case_count_final = 0
                for i in range(len(case_dic)):
                    case_count_final += len(case_dic[f'{i + 1}'][1])

                dic_count[f'sg_{sym_no}'] = case_count_final
            elif atom_type_no == 2:
                case_count_final = 0
                
                for n in range(1, len(case_dic) + 1):
                    dic_check_list = case_dic[f'{n}'][1]
                    dic_remove_list = []
                    for m in range(len(case_dic[f'{n}'][1])):
                        wyckoff_list = [case_dic[f'{n}'][1][m]]
                        check_list = []
                        for l in range(len(case_dic[f'{n}'][1][m])):
                            wyckoff_symble = case_dic[f'{n}'][1][m][l]
                            wyckoff_muti = sg[f'sg_{sym_no}'][wyckoff_symble][0]
                            check_add = [wyckoff_symble, wyckoff_muti]
                            check_list.append(check_add)
                            # print(f'check_add:{check_add}')
                        case_count_inter = divide_list_to_equal_sum_count(check_list, self.ratio)
                        case_count_final += case_count_inter
                
                dic_count[f'sg_{sym_no}'] = case_count_final                

        return dic_count

    def case_list(self, sym_no, mode):
        case_dic = {}
        case_no = 0
        # single atom type number
        atom_type_no = len(self.ratio) - 1
        sg_len = len(sg[f'sg_{sym_no}'])
        # generate muti_list
        muti_list = []
        for m in range(1, sg_len):
            muti_list.append([f's{m}', int(sg[f'sg_{sym_no}'][f's{m}'][0])])
            # muti_list.append(int(sg[f'sg_{sym_no}'][f's{m}'][0]))

        # rigid_at_general = 0
        case_count = 0
        sg_range = 0
        if self.gen == 0:
            sg_range = sg_len
        elif self.gen == 1:
            sg_range = sg_len + 1
        # print(f'sg_range:{sg_range}')
        wyckoff_sym = rigid_wyckoff_sym(self.rigid_type, sym_no)
        for i in range(1, sg_range):
            if sg[f'sg_{sym_no}'][f's{i}'][1] in wyckoff_sym:
                site_muti = int(sg[f'sg_{sym_no}'][f's{i}'][0])
                if site_muti <= self.rigid_no:
                    ratio_total = 0
                    for j in range(1, len(self.ratio)):
                        ratio_total += self.ratio[j]
                    single_atom_no = site_muti * ratio_total
                    # the max number of atoms can be at general position
                    gen_max = int(single_atom_no / int(sg[f'sg_{sym_no}'][f's{sg_len}'][0]))
                    gen_muti = sg[f'sg_{sym_no}'][f's{sg_len}'][0]

                    if sg[f'sg_{sym_no}'][f's{i}'][3] == 1:
                        # print(f'muti_list:{muti_list}')
                        for gen_no in range(gen_max + 1):
                            # for n in range(0, len(muti_list) + 1):
                            for n in range(0, single_atom_no):
                                count2 = 0
                                enum = list(combinations_with_replacement(muti_list, n))
                                enum_inter = enum.copy()
                                # print(enum_inter)
                                if mode == 0 and len(enum) > 100000:
                                    case_count = 'TBD'
                                    break
                                elif mode == 1 or (mode == 0 and len(enum) <= 100000):
                                    if len(enum) != 0:
                                        len_k = len(enum[0])
                                        for k in enum:
                                            # print(f'k:{k}')
                                            sum_k = 0
                                            for k1 in range(len_k):
                                                sum_k += int(k[k1][1])
                                            if sum_k + (gen_muti * gen_no) != single_atom_no:
                                                enum_inter.remove(k)
                                    # print(f'rigid_position:{i}')
                                    # print(f'atom_at_general:{gen_no}')
                                    # print(f'n:{n}')
                                    # print(enum)
                                    # print(enum_inter)
                                    # print(f'sum:{sum_k + (gen_muti * gen_no)}, single_atom_no:{single_atom_no}')
                                    # print('-----------')
                                    if len(enum_inter) != 0:
                                        case_no += 1
                                        case_list = []
                                        for k2 in range(len(enum_inter)):
                                            case_list_inter = []
                                            for k3 in range(len_k):
                                                case_list_inter.append(enum_inter[k2][k3][0])
                                            for k4 in range(gen_no):
                                                case_list_inter.append(f's{sg_len}')
                                            case_list.append(case_list_inter)
                                        case_dic[f'{case_no}'] = [f'{i}', case_list]

                                        count2 += 1
                                        # print(f'sg:{sym_no}')
                                        # print(f'site:s{i}')
                                        # print(f'gen number:{gen_no}')
                                        # print(f'atom at special:{n}')
                                        # print(f'expected:{len(enum)}')
                                        # print(f'current:{count2}')
                                        # print('-------------------')
                                    enum_len = len(enum_inter)
                                    if atom_type_no == 1:
                                        case_count += enum_len
                                        # print(f'case count:{case_count}')
                                    elif atom_type_no == 2:
                                        for k in range(enum_len):
                                            enum_list = list(enum[k])
                                            check_atom = divide_list_to_equal_sum(enum_list, self.ratio)
                                            if check_atom == 'yes':
                                                case_count += 1
                            if case_count == 'TBD':
                                break
                        if case_count == 'TBD':
                            break
        # get rid of atom at fixed wyckoff position
        '''
        (list become empty)
        '''
        s_list = []
        for i in range(sg_len):
            s_x = sg[f'sg_{sym_no}'][f's{i + 1}'][2][0]
            s_y = sg[f'sg_{sym_no}'][f's{i + 1}'][2][1]
            s_z = sg[f'sg_{sym_no}'][f's{i + 1}'][2][2]
            if s_x != 'x' and s_y != 'y' and s_z != 'z':
                s_list.append(f's{i + 1}')

        s_list_len = len(s_list)
        # print(s_list)
        if s_list_len != 0:
            dic_len = len(case_dic)
            for i in range(dic_len):
                check_list = case_dic[f'{i + 1}'][1]
                # print(f'check_list:{check_list}')
                # print(f'check_list_len:{len(check_list)}')
                j_remove = []
                for j in range(len(check_list)):
                    s_check = 0
                    for k in range(s_list_len):
                        if check_list[j].count(s_list[k]) >= 2:
                            # print(f'j:{j}')
                            # print(f'{s_list[k]} count:{check_list[j].count(s_list[k])}')
                            s_check = 1
                            # print(f's_check:{s_check}')
                            # print('------------------------------')
                            break
                    if s_check == 1:
                        j_remove.append(j)
                j_remove_len = len(j_remove)
                print(f'j_remove:{j_remove}')
                j_remove_list = []
                if j_remove_len != 0:
                    for n in j_remove:
                        j_remove_list.append(check_list[n])
                for m in j_remove_list:
                    check_list.remove(m)
                # print(case_dic[f'{i + 1}'][1])
                # print(check_list)
                # case_dic[f'{i}'][1] = check_list
        return case_dic


class SymCasesP1:
    """
    ---------------------------------------------------------------
    todo list:
    mutiple rigid body before sym transform
    mutiple type of rigid bodies
    combine case_count and case_list
    ---------------------------------------------------------------
    ratio = array
    rigid_at_general = 1: rigid body allow to be at general position
                     = 0: not allowd to be at general position
    rigid_number: the max number of rigid bodies
    eg: sym = SymCases(tetrahedron, np.array([1, 3]), 1, 3)
    """    
    def __init__(self, rigid_type, ratio, rigid_at_general, rigid_number):
        self.rigid_type = rigid_type
        self.ratio = ratio
        self.gen = rigid_at_general
        self.rigid_no = rigid_number

    """
    mode = 1: count all number no matter how large
         = 0: TBD if number is large than 100000
    """
    def case_list(self):
        case_dic = {}
        for i in range(self.rigid_no):
            rigid_list = []
            single_list = []
            for j in range(i + 1):
                rigid_list.append(1)
            for k in range((i + 1) * 3):
                single_list.append('s1')
            case_dic[f'{i + 1}'] = [rigid_list, single_list]
        return case_dic
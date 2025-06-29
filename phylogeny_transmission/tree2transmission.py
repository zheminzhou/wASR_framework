import ete3_extensions, click


@click.command()
@click.option('-p', '--proportion', help='proportion for a state to form transmission', show_default=True, default=0.5, type=float)
@click.option('-a', '--all_nodes', help='include information for nodes without transmissions', default=False, show_default=True, is_flag=True)
@click.argument('nexus')
def main(nexus, proportion, all_nodes) :
    trees = ete3_extensions.read_nexus(nexus)
    for tre in trees :
        for node in tre.get_descendants('postorder') :
            if node.is_leaf() :
                node.d = 1
            else :
                node.d = sum([c.d for c in node.children ])

            parent = node.up
            if all_nodes :
                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(parent.annotations['state'], node.annotations['state'], parent.annotations['state.prop'], node.annotations['state.prop'], \
                    parent.annotations['date'], node.annotations['date'], node.d, parent.name, node.name))
            else :
                if node.annotations['state.prop'] < proportion :
                    continue
                while parent and parent.annotations['state.prop'] < proportion and parent.up :
                    if parent.annotations['state'] != parent.up.annotations['state'] and parent.annotations['state'] != node.annotations['state'] :
                        break
                    parent = parent.up
                
                if parent and parent.annotations['state'] != node.annotations['state'] and parent.annotations['state'] != 'ND' and node.annotations['state'] != 'ND' :
                    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(parent.annotations['state'], node.annotations['state'], parent.annotations['state.prop'], node.annotations['state.prop'], \
                        parent.annotations['date'], node.annotations['date'], node.d, parent.name, node.name))
    




if __name__ == '__main__' :
    main()

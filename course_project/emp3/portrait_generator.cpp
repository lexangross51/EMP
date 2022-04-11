#include "portrait_generator.hpp"

void portrait_generator::portrait(space_grid& grid, std::vector<uint32_t>& ig, std::vector<uint32_t>& jg)
{
    std::vector<std::set<uint32_t>> list(grid.get_nodes_count());

    for (uint32_t ielem = 0; ielem < grid.get_elems_count(); ielem++)
    {
        finite_element elem = grid.get_elem(ielem);

        for (uint32_t i = 0; i < 8; i++)
        {
            uint32_t elem_to_insert = elem[i];

            for (uint32_t j = i + 1; j < 9; j++)
            {
                uint32_t insert_pos = elem[j];

                list[insert_pos].insert(elem_to_insert);
            }
        }
    }
    
    ig.resize(grid.get_nodes_count() + 1);

    ig[0] = 0;
    ig[1] = 0;
    for (size_t i = 1; i < list.size(); i++)
        ig[i + 1] = ig[i] + list[i].size();

    jg.resize(ig.back());

    for (size_t i = 1, j = 0; i < list.size(); i++)
    {
        for (const auto& it : list[i])
            jg[j++] = it;
    }
}
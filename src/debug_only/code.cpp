Put node_start in the OPEN list with f (node_start) = h(node_start) (initialization)

while the OPEN list is not empty {
Take from the open list the node node_current with the lowest
f (node_current) = g(node_current) + h(node_current)

if node_current is node_goal we have found the solution; break

Generate each state node_successor that come after node_current
for each node_successor of node_current 
{
		Set successor_current_cost = g(node_current) + w(node_current, node_successor)
		
		if node_successor is in the OPEN list 
		{
			if g(node_successor) ≤ successor_current_cost continue (to line 31 )
		} 
		else if node_successor is in the CLOSED list 
		{
			if g(node_successor) ≤ successor_current_cost continue (to line 31 )
			Move node_successor from the CLOSED list to the OPEN list
		} 
		else 
		{
			Add node_successor to the OPEN list
			Set h(node_successor) to be the heuristic distance to node_goal
		}
		
		Set g(node_successor) = successor_current_cost
		Set the parent of node_successor to node_current
}

Add node_current to the CLOSED list
}

if(node_current != node_goal) exit with error (the OPEN list is empty)
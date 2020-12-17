/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#pragma once

#include <memory>
#include <optional>
#include <string>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../traits.hpp"
#include "../utils/algorithm.hpp"
#include "detail/foreach.hpp"

namespace mockturtle
{

class abstract_mig_network
{
public:
#pragma region Types and constructors
  static constexpr auto min_fanin_size = 3u; 
  static constexpr auto max_fanin_size = std::numeric_limits<uint32_t>::max();

  using base_type = abstract_mig_network;
  using node = uint32_t;

/*
  A signal represents an edge in the network
*/
  struct signal
  {
    signal() = default;
    signal(uint32_t index, bool complement=false):index(index), complement(complement)
    {}

    uint32_t index;   // index of the node that this edge is pointing to
    bool complement;  // denoting whether the edge is complemented or not

    signal operator!() const
    { 
      return {index, !complement};
    }

    signal operator+() const
    { 
      return {index, false};  // return the non-complemented signal
    }

    signal operator-() const
    { 
      return {index, true}; // return the complemented signal
    }

    signal operator^( bool complement ) const
    {  
      return {index, this->complement != complement}; // return the signal inverting the complemented
    }

    bool operator==( signal const& other ) const
    {  
      return index == other.index && complement == other.complement;
    }

    bool operator!=( signal const& other ) const
    {  
      return index != other.index || complement != other.complement;
    }

    bool operator<( signal const& other ) const
    {  
      return index < other.index || (index == other.index && !complement && other.complement);
    }
    
  };

  struct storage_type
  {
    struct node_type {
      std::vector<signal> fanin{}; // saves the input signals to a particular node  
      uint32_t fanout_size{};  // number of outputs of the node
    };   // represents a single node in the graph

    storage_type()
    {
      /* constant 0 node */
      nodes.emplace_back(); 
    }

    std::vector<node_type> nodes;
    std::vector<uint32_t> inputs;  // indices of input nodes 
    std::vector<signal> outputs; // saves all output signals of the network
  };
  using storage = std::shared_ptr<storage_type>;

  abstract_mig_network()
    : _storage( std::make_shared<storage_type>() )
  {}

  abstract_mig_network(std::shared_ptr<storage_type> storage)
    :_storage(storage)
  {}

#pragma endregion

#pragma region Primary I / O and constants
  signal get_constant( bool value ) const
  {
    return signal{0, value};  // index of a constant is always 0
  }

  signal create_pi( std::string const& name = std::string() )
  {
    (void)name;  // ignore the string name for now

    // calculate the index from the size of the nodes vector
    const auto index = static_cast<uint32_t>( _storage->nodes.size() );

    // add a new empty entry to nodes vector, because it has no fanin and fanout is still unknown at that moment
    _storage->nodes.emplace_back(); 

    // the node is an input so add its index to the nodes vector
    _storage->inputs.emplace_back( index );  // Since inputs vector holds the indices of input nodes

    // return the signal (edge) coming out of your newly created node
    return {index, 0}; // primary inputs aren't complemented initially
  }

  uint32_t create_po( signal const& f, std::string const& name = std::string() )
  {
    (void)name;  // ignore the string name for now

    /* f denotes the output signal, the index in f denotes the node that this output will emerge from 
       therefore, increase the fanout of the node at that index
    */
    _storage->nodes[f.index].fanout_size++;

    // create an entry for this new output signal in outputs vector
    auto const po_index = static_cast<uint32_t>( _storage->outputs.size() );
    _storage->outputs.emplace_back( f );

    return po_index;
  }

  bool is_pi( node const& n ) const
  {
    // note that node n denotes just an index of the node under question

    // n has to be > 0  because index 0 is preserved for the constant
    return n > 0 && _storage->nodes[n].fanin.size() == 0u;
  }

  bool constant_value( node const& n ) const
  {
    (void)n;
    return false;     
  }
#pragma endregion

#pragma region Create binary functions
  signal create_and( signal a, signal b )
  {
    // since this maj will be only true if both a and b are true which is AND
    return create_maj( get_constant( false ), a, b );
  }

  signal create_or( signal const& a, signal const& b )
  {
    // TODO
    // since this maj will be true if at least a or b are true which is OR
    return create_maj( get_constant( true ), a, b );
  }

  signal create_xor( signal const& a, signal const& b )
  {
    const auto fcompl = a.complement ^ b.complement;
    const auto c1 = create_and( +a, -b );
    const auto c2 = create_and( +b, -a );
    return create_and( !c1, !c2 ) ^ !fcompl;
  }
#pragma endregion

#pragma region Create ternary functions
  signal create_maj( signal const& a, signal const& b, signal const& c )
  {
    std::vector<signal> fs;
    fs.emplace_back(a);
    fs.emplace_back(b);
    fs.emplace_back(c);

    return create_nary_maj( fs );
  }
#pragma endregion

#pragma region Create nary functions
  signal create_nary_maj( std::vector<signal> const& fs )
  {
    /*
      if number of signals isn't odd, raise an error
      if condition is false (when fs.size() is odd) an error is raised
    */
    assert((fs.size() % 2)==1);

    const auto index = static_cast<uint32_t>( _storage->nodes.size() );

    storage_type::node_type new_gate;
    for(int i = 0; i < fs.size(); i++)
    {
      new_gate.fanin.emplace_back(fs.at(i));
    }

    _storage->nodes.emplace_back(new_gate);

    // return the signal (edge) coming out of your newly created maj gate node
    return {index, 0}; 

  }
#pragma endregion

#pragma region Nodes and signals
  node get_node( signal const& f ) const
  {
    // return the index stored in the passed signal f
    return f.index;
  }

  bool is_complemented( signal const& f ) const
  {
    return f.complement;
  }

  uint32_t node_to_index( node const& n ) const
  {
    return static_cast<uint32_t>( n );
  }
#pragma endregion

#pragma region Node and signal iterators
  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( _storage->inputs.begin(), _storage->inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    detail::foreach_element( _storage->outputs.begin(), _storage->outputs.end(), fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint32_t>( 1u, _storage->nodes.size() ); /* start from 1 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_pi( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_fanin( node const& n, Fn&& fn ) const
  {
    if ( n == 0 || is_pi( n ) )
      return;

    // iterate over the fanin nodes of node n
    detail::foreach_element( _storage->nodes.at(n).fanin.begin(), _storage->nodes.at(n).fanin.end(), fn );
  }
#pragma endregion

#pragma region Structural properties
  uint32_t size() const
  {
    return _storage->nodes.size(); // total number of nodes 
  }

  uint32_t num_pis() const
  {
    return static_cast<uint32_t>( _storage->inputs.size() );  // total number of inputs
  }

  uint32_t num_pos() const
  {
    return static_cast<uint32_t>( _storage->outputs.size() );  // total number of outputs
  }

  uint32_t num_gates() const
  {
    // size includes the maj gates + inputs + constant value that's why subtract the extra 1
    return size() - num_pis() - 1;  // total number of majority gates  
  }

  uint32_t fanin_size( node const& n ) const
  {
    return _storage->nodes[n].fanin.size();
  }

  uint32_t fanout_size( node const& n ) const
  {
    return _storage->nodes[n].fanout_size;
  }

#pragma endregion

#pragma region Value simulation
  template<typename Iterator>
  iterates_over_truth_table_t<Iterator>
  compute( node const& n, Iterator begin, Iterator end ) const
  {
    (void) end;
    assert( n != 0 && !is_pi( n ) );
    typename Iterator::value_type maj_n(fanin_size(n));
    kitty::create_majority(maj_n);
    std::vector<typename Iterator::value_type> tts;
    foreach_fanin(n, [&](signal s, uint32_t i) {
      tts.push_back( is_complemented(s) ? ~(*(begin + i)) : *(begin + i));
    });
    return kitty::compose_truth_table(maj_n, tts);
  }
#pragma endregion

private:
  storage _storage;  // this is a shared_ptr pointing to storage_type struct

};

} // namespace mockturtle

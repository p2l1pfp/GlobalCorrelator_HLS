library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity tk_router_element is
    port(
            ap_clk : IN STD_LOGIC;
            enabled : IN STD_LOGIC;
            newevent : IN STD_LOGIC;
            links_in : IN particles(NTKFIBERS-1 downto 0);
            fifo_same : OUT particles(NTKFIBERS-1 downto 0);
            fifo_next : OUT particles(NTKFIBERS-1 downto 0);
            fifo_prev : OUT particles(NTKFIBERS-1 downto 0);
            fifo_same_write: OUT std_logic_vector(NTKFIBERS-1 downto 0);
            fifo_next_write: OUT std_logic_vector(NTKFIBERS-1 downto 0);
            fifo_prev_write: OUT std_logic_vector(NTKFIBERS-1 downto 0);
            fifo_same_roll: OUT std_logic_vector(NTKFIBERS-1 downto 0);
            fifo_next_roll: OUT std_logic_vector(NTKFIBERS-1 downto 0);
            fifo_prev_roll: OUT std_logic_vector(NTKFIBERS-1 downto 0)
    );
end tk_router_element;

architecture Behavioral of tk_router_element is
begin
    link2fifo : process(ap_clk)
        variable link_this, link_next, link_prev : std_logic;
    begin
        if rising_edge(ap_clk) then
            for ifib in 0 to NTKFIBERS-1 loop
                if enabled = '0' or links_in(ifib).pt = 0 then
                    link_this := '0';
                    link_prev := '0';
                    link_next := '0';
                else
                    link_this := '1';
                    if links_in(ifib).phi >= PHI_MARGIN_POS then
                        link_prev := '0';
                        link_next := '1';
                    elsif links_in(ifib).phi <= PHI_MARGIN_NEG then
                        link_prev := '1';
                        link_next := '0';
                    else
                        link_prev := '0';
                        link_next := '0';
                    end if;
                end if;
                fifo_same(ifib)      <= links_in(ifib);
                fifo_next(ifib).pt   <= links_in(ifib).pt;
                fifo_next(ifib).eta  <= links_in(ifib).eta;
                fifo_next(ifib).phi  <= links_in(ifib).phi - PHI_SHIFT;
                fifo_next(ifib).rest <= links_in(ifib).rest;
                fifo_prev(ifib).pt   <= links_in(ifib).pt;
                fifo_prev(ifib).eta  <= links_in(ifib).eta;
                fifo_prev(ifib).phi  <= links_in(ifib).phi + PHI_SHIFT;
                fifo_prev(ifib).rest <= links_in(ifib).rest;
                fifo_same_write(ifib) <= link_this;
                fifo_next_write(ifib) <= link_next;
                fifo_prev_write(ifib) <= link_prev;
                fifo_same_roll(ifib)  <= newevent;
                fifo_next_roll(ifib)  <= newevent;
                fifo_prev_roll(ifib)  <= newevent;
            end loop;
        end if;
    end process link2fifo;

end Behavioral;
